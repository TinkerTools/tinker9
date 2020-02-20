#include "elec.h"
#include "io_fort_str.h"
#include "md.h"
#include "pme.h"
#include "potent.h"
#include <tinker/detail/chgpot.hh>
#include <tinker/detail/limits.hh>
#include <tinker/detail/mpole.hh>

TINKER_NAMESPACE_BEGIN
int use_elec()
{
   return use_potent(mpole_term) || use_potent(polar_term);
}

int use_ewald()
{
   return limits::use_ewald;
}

static void pole_data_(rc_op op)
{
   if (op & rc_dealloc) {
      device_array::deallocate(zaxis, pole, rpole, udir, udirp, uind, uinp);

      device_array::deallocate(trqx, trqy, trqz, vir_trq);
   }

   if (op & rc_alloc) {
      device_array::allocate(n, &zaxis, &pole, &rpole);

      if (use_potent(polar_term)) {
         device_array::allocate(n, &uind, &uinp, &udir, &udirp);
      } else {
         uind = nullptr;
         uinp = nullptr;
         udir = nullptr;
         udirp = nullptr;
      }

      if (rc_flag & calc::grad) {
         device_array::allocate(n, &trqx, &trqy, &trqz);
      } else {
         trqx = nullptr;
         trqy = nullptr;
         trqz = nullptr;
      }

      if (rc_flag & calc::virial) {
         device_array::allocate(buffer_size(), &vir_trq);
         virial_buffers.push_back(vir_trq);
      } else {
         vir_trq = nullptr;
      }
   }

   if (op & rc_init) {
      electric = chgpot::electric;
      dielec = chgpot::dielec;

      // Regarding chkpole routine:
      // 1. The chiralities of the atoms will not change in the simulations;
      // 2. chkpole routine has been called in mechanic routine so that the
      // values in mpole::pole are correct;
      // 3. yaxis values are directly copied from Tinker, and are NOT
      // subtracted by 1 becasue of the checks in chkpole;
      // 4. GPU chkpole kernel is necessary when unexpected changes of
      // charalities may happen, e.g. in Monte Carlo simulations.
      static_assert(sizeof(LocalFrame) == 4 * sizeof(int), "");
      std::vector<LocalFrame> zaxisbuf(n);
      for (int i = 0; i < n; ++i) {
         zaxisbuf[i].zaxis = mpole::zaxis[i] - 1;
         zaxisbuf[i].xaxis = mpole::xaxis[i] - 1;
         zaxisbuf[i].yaxis = mpole::yaxis[i];
         fstr_view str = mpole::polaxe[i];
         int val;
         if (str == "Z-Only")
            val = pole_z_only;
         else if (str == "Z-then-X")
            val = pole_z_then_x;
         else if (str == "Bisector")
            val = pole_bisector;
         else if (str == "Z-Bisect")
            val = pole_z_bisect;
         else if (str == "3-Fold")
            val = pole_3_fold;
         else
            val = pole_none;
         zaxisbuf[i].polaxe = val;
      }
      device_array::copyin(WAIT_NEW_Q, n, zaxis, zaxisbuf.data());

      std::vector<double> polebuf(mpl_total * n);
      for (int i = 0; i < n; ++i) {
         int b1 = mpl_total * i;
         int b2 = mpole::maxpole * i;
         // Tinker c = 0, dx = 1, dy = 2, dz = 3
         // Tinker qxx = 4, qxy = 5, qxz = 6
         //        qyx    , qyy = 8, qyz = 9
         //        qzx    , qzy    , qzz = 12
         polebuf[b1 + mpl_pme_0] = mpole::pole[b2 + 0];
         polebuf[b1 + mpl_pme_x] = mpole::pole[b2 + 1];
         polebuf[b1 + mpl_pme_y] = mpole::pole[b2 + 2];
         polebuf[b1 + mpl_pme_z] = mpole::pole[b2 + 3];
         polebuf[b1 + mpl_pme_xx] = mpole::pole[b2 + 4];
         polebuf[b1 + mpl_pme_xy] = mpole::pole[b2 + 5];
         polebuf[b1 + mpl_pme_xz] = mpole::pole[b2 + 6];
         polebuf[b1 + mpl_pme_yy] = mpole::pole[b2 + 8];
         polebuf[b1 + mpl_pme_yz] = mpole::pole[b2 + 9];
         polebuf[b1 + mpl_pme_zz] = mpole::pole[b2 + 12];
      }
      device_array::copyin(WAIT_NEW_Q, n, pole, polebuf.data());
   }
}

void elec_data(rc_op op)
{
   if (!use_elec())
      return;

   rc_man pole42_{pole_data_, op};
   rc_man pme42_{pme_data, op};
}

extern void chkpole();
extern void rotpole();
void elec_init(int vers)
{
   if (!use_elec())
      return;

   if (vers & calc::grad)
      device_array::zero(PROCEED_NEW_Q, n, trqx, trqy, trqz);
   if (vers & calc::virial)
      device_array::zero(PROCEED_NEW_Q, buffer_size(), vir_trq);

   chkpole();
   rotpole();

   if (use_ewald()) {
      pme_init(vers);
#if TINKER_CUDART
      if (epme_unit.valid()) {
         bspline_fill(epme_unit, 3);
      }


      if (ppme_unit.valid() && (ppme_unit != epme_unit)) {
         bspline_fill(ppme_unit, 2);
      }


      if (pvpme_unit.valid()) {
         bspline_fill(pvpme_unit, 2);
      }
#endif
   }
}
TINKER_NAMESPACE_END
