#include "elec.h"

#include "ext/tinker/detail/chgpot.hh"
#include "ext/tinker/detail/limits.hh"
#include "ext/tinker/detail/mpole.hh"
#include "io_fort_str.h"
#include "md.h"
#include "pme.h"
#include "potent.h"

TINKER_NAMESPACE_BEGIN
int use_elec() { return use_potent(mpole_term) || use_potent(polar_term); }

int use_ewald() { return limits::use_ewald; }

static void pole_data_(rc_op op) {
  if (op & rc_dealloc) {
    DeviceMemory::deallocate_bytes(zaxis);

    vir_trq_handle.close();
  }

  if (op & rc_alloc) {
    DeviceMemory::allocate_bytes(&zaxis, n * sizeof(LocalFrame));
    pole_vec.reserve(n * mpl_total);
    rpole_vec.reserve(n * mpl_total);

    if (use_potent(polar_term)) {
      uind_vec.reserve(3 * n);
      uinp_vec.reserve(3 * n);
      udir_vec.reserve(3 * n);
      udirp_vec.reserve(3 * n);
    } else {
      uind_vec.clear();
      uinp_vec.clear();
      udir_vec.clear();
      udirp_vec.clear();
    }

    if (rc_flag & calc::grad) {
      trqx_vec.reserve(n);
      trqy_vec.reserve(n);
      trqz_vec.reserve(n);
    } else {
      trqx_vec.clear();
      trqy_vec.clear();
      trqz_vec.clear();
    }

    vir_trq_handle = Virial::open();
    vir_trq_handle->alloc(n);
  }

  if (op & rc_init) {
    electric = chgpot::electric;
    dielec = chgpot::dielec;

    // Regarding chkpole routine:
    // 1. The chiralities of the atoms will not change in the simulations;
    // 2. chkpole routine has been called in mechanic routine so that the values
    // in mpole::pole are correct;
    // 3. yaxis values are directly copied from Tinker, and are NOT
    // subtracted by 1 becasue of the checks in chkpole;
    // 4. GPU chkpole kernel is necessary when unexpected changes of charalities
    // may happen, e.g. in Monte Carlo simulations.
    static_assert(sizeof(LocalFrame) == 4 * sizeof(int), "");
    std::vector<int> zaxisbuf(4 * n);
    for (int i = 0; i < n; ++i) {
      int base = 4 * i;
      zaxisbuf[base] = mpole::zaxis[i] - 1;
      zaxisbuf[base + 1] = mpole::xaxis[i] - 1;
      zaxisbuf[base + 2] = mpole::yaxis[i];
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
      zaxisbuf[base + 3] = val;
    }
    DeviceMemory::copyin_array(reinterpret_cast<int*>(zaxis), zaxisbuf.data(),
                               4 * n);

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
    pole_vec.copyin(polebuf.data(), mpl_total * n);
  }
}

void elec_data(rc_op op) {
  if (!use_elec())
    return;

  rc_man pole42_{pole_data_, op};
  rc_man pme42_{pme_data, op};
}

extern void chkpole();
extern void rotpole();
void elec_init(int vers) {
  if (!use_elec())
    return;

  // zero torque

  if (vers & calc::grad) {
    trqx_vec.zero(n);
    trqy_vec.zero(n);
    trqz_vec.zero(n);
  }

  // zero torque-related virial

  if (vers & calc::virial) {
    vir_trq_handle->zero();
  }

  chkpole();
  rotpole();

  if (use_ewald())
    pme_init(vers);
}
TINKER_NAMESPACE_END
