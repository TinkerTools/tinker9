#include "elec.h"
#include "io_fort_str.h"
#include "md.h"
#include "pmestuf.h"
#include "potent.h"
#include <tinker/detail/chgpot.hh>
#include <tinker/detail/limits.hh>
#include <tinker/detail/mpole.hh>


TINKER_NAMESPACE_BEGIN
real electric, dielec;


bool use_ewald()
{
   return limits::use_ewald;
}


//====================================================================//


real* pchg;


void pchg_data(rc_op op)
{
   if (!use_potent(charge_term))
      return;
}


//====================================================================//


LocalFrame* zaxis;
real (*pole)[mpl_total];
real (*rpole)[mpl_total];
real *trqx, *trqy, *trqz;
virial_buffer vir_trq;
real (*udir)[3];
real (*udirp)[3];
real (*uind)[3];
real (*uinp)[3];


void pole_data(rc_op op)
{
   if (!use_potent(mpole_term) && !use_potent(polar_term))
      return;


   if (op & rc_dealloc) {
      darray::deallocate(zaxis, pole, rpole, udir, udirp, uind, uinp);
      darray::deallocate(trqx, trqy, trqz, vir_trq);
   }


   if (op & rc_alloc) {
      darray::allocate(n, &zaxis, &pole, &rpole);


      if (use_potent(polar_term)) {
         darray::allocate(n, &uind, &uinp, &udir, &udirp);
      } else {
         uind = nullptr;
         uinp = nullptr;
         udir = nullptr;
         udirp = nullptr;
      }


      if (rc_flag & calc::grad) {
         darray::allocate(n, &trqx, &trqy, &trqz);
      } else {
         trqx = nullptr;
         trqy = nullptr;
         trqz = nullptr;
      }
      if (rc_flag & calc::virial) {
         darray::allocate(buffer_size(), &vir_trq);
         virial_buffers.push_back(vir_trq);
      } else {
         vir_trq = nullptr;
      }
   }


   if (op & rc_init) {
      // Regarding chkpole routine:
      // 1. The chiralities of the atoms are unlikely to change in MD; but
      // still possible in Monte Carlo;
      // 2. yaxis values are directly copied from Tinker, and are NOT
      // subtracted by 1 becasue of the checks in chkpole.
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
      darray::copyin(WAIT_NEW_Q, n, zaxis, zaxisbuf.data());


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
      darray::copyin(WAIT_NEW_Q, n, pole, polebuf.data());
   }
}


//====================================================================//


void elec_data(rc_op op)
{
   if (op & rc_init) {
      electric = chgpot::electric;
      dielec = chgpot::dielec;
   }
   rc_man pchg42{pchg_data, op};
   rc_man pole42{pole_data, op};
   rc_man pme42{pme_data, op};
}


//====================================================================//


void mpole_init(int vers)
{
   if (!use_potent(mpole_term) && !use_potent(polar_term))
      return;


   if (vers & calc::grad)
      darray::zero(PROCEED_NEW_Q, n, trqx, trqy, trqz);
   if (vers & calc::virial)
      darray::zero(PROCEED_NEW_Q, buffer_size(), vir_trq);


   chkpole();
   rotpole();


   if (use_ewald()) {
      rpole_to_cmp();
      if (vir_m)
         darray::zero(PROCEED_NEW_Q, buffer_size(), vir_m);
      if (pltfm_config & CU_PLTFM) {
         bool precompute_theta = (!TINKER_CU_THETA_ON_THE_FLY_GRID_MPOLE) ||
            (!TINKER_CU_THETA_ON_THE_FLY_GRID_UIND);
         if (epme_unit.valid()) {
            if (precompute_theta)
               bspline_fill(epme_unit, 3);
         }
         if (ppme_unit.valid() && (ppme_unit != epme_unit)) {
            if (precompute_theta)
               bspline_fill(ppme_unit, 2);
         }
         if (pvpme_unit.valid()) {
            if (precompute_theta)
               bspline_fill(pvpme_unit, 2);
         }
      }
   }
}


void chkpole()
{
   chkpole_acc();
}


void rotpole()
{
   rotpole_acc();
}


void torque(int vers)
{
   // #if TINKER_CUDART
   //    if (pltfm_config & CU_PLTFM)
   //       torque_cu(vers);
   //    else
   // #endif
   torque_acc(vers);
}
TINKER_NAMESPACE_END
