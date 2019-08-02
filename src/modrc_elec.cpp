#include "gpu/e_mpole.h"
#include "gpu/e_polar.h"
#include "mod_elec.h"
#include "mod_md.h"
#include "mod_pme.h"
#include "util_array.h"
#include "util_potent.h"
#include "util_rt.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
int use_elec() { return use_potent(mpole_term) || use_potent(polar_term); }

void elec_data(rc_op op) {
  if (!use_elec())
    return;

  if (op & rc_dealloc) {
    dealloc_array(zaxis);
    dealloc_array(pole);
    dealloc_array(rpole);

    dealloc_array(uind);
    dealloc_array(uinp);
    dealloc_array(udir);
    dealloc_array(udirp);

    dealloc_array(trqx);
    dealloc_array(trqy);
    dealloc_array(trqz);
    dealloc_array(vir_trq);
  }

  if (op & rc_alloc) {
    const size_t rs = sizeof(real);
    size_t size;

    size = sizeof(local_frame_t);
    alloc_array(&zaxis, n * size);
    size = rs * mpl_total;
    alloc_array(&pole, n * size);
    alloc_array(&rpole, n * size);

    if (use_potent(polar_term)) {
      alloc_array(&uind, 3 * n * rs);
      alloc_array(&uinp, 3 * n * rs);
      alloc_array(&udir, 3 * n * rs);
      alloc_array(&udirp, 3 * n * rs);
    } else {
      uind = nullptr;
      uinp = nullptr;
      udir = nullptr;
      udirp = nullptr;
    }

    if (use_data & calc::grad) {
      alloc_array(&trqx, rs * n);
      alloc_array(&trqy, rs * n);
      alloc_array(&trqz, rs * n);
    } else {
      trqx = nullptr;
      trqy = nullptr;
      trqz = nullptr;
    }

    alloc_array(&vir_trq, rs * 9);
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
    static_assert(sizeof(local_frame_t) == 4 * sizeof(int), "");
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
    copyin_array(reinterpret_cast<int*>(zaxis), zaxisbuf.data(), 4 * n);

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
    copyin_array(reinterpret_cast<real*>(pole), polebuf.data(), mpl_total * n);
  }

  pme_data(op);
}

extern void chkpole();
extern void rotpole();
void elec_init(int vers) {
  if (!use_elec())
    return;

  // zero torque

  if (vers & calc::grad) {
    zero_array(trqx, n);
    zero_array(trqy, n);
    zero_array(trqz, n);
  }

  // zero torque-related virial

  if (vers & calc::virial) {
    zero_array(vir_trq, 9);
  }

  chkpole();
  rotpole();

  if (use_ewald()) {
    pme_init(vers);
  }
}
TINKER_NAMESPACE_END
