#include "gpu/decl_elec.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_pme.h"
#include "gpu/e_mpole.h"
#include "gpu/e_polar.h"
#include "gpu/rc.h"
#include "util/fort_str.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
double mpole_switch_cut, mpole_switch_off;

local_frame_t* zaxis;

real (*pole)[mpl_total];
real (*rpole)[mpl_total];

real (*uind)[3];
real (*uinp)[3];
real (*udir)[3];
real (*udirp)[3];

real *trqx, *trqy, *trqz;
real* vir_trq;

int use_elec() { return use_empole() || use_epolar(); }

void elec_data(rc_t rc) {
  if (!use_elec())
    return;

  if (rc & rc_dealloc) {
    check_cudart(cudaFree(zaxis));
    check_cudart(cudaFree(pole));
    check_cudart(cudaFree(rpole));

    check_cudart(cudaFree(uind));
    check_cudart(cudaFree(uinp));
    check_cudart(cudaFree(udir));
    check_cudart(cudaFree(udirp));

    check_cudart(cudaFree(trqx));
    check_cudart(cudaFree(trqy));
    check_cudart(cudaFree(trqz));
    check_cudart(cudaFree(vir_trq));
  }

  if (rc & rc_alloc) {
    const size_t rs = sizeof(real);
    size_t size;

    size = sizeof(local_frame_t);
    check_cudart(cudaMalloc(&zaxis, n * size));
    size = rs * mpl_total;
    check_cudart(cudaMalloc(&pole, n * size));
    check_cudart(cudaMalloc(&rpole, n * size));

    if (use_epolar()) {
      check_cudart(cudaMalloc(&uind, 3 * n * rs));
      check_cudart(cudaMalloc(&uinp, 3 * n * rs));
      check_cudart(cudaMalloc(&udir, 3 * n * rs));
      check_cudart(cudaMalloc(&udirp, 3 * n * rs));
    } else {
      uind = nullptr;
      uinp = nullptr;
      udir = nullptr;
      udirp = nullptr;
    }

    if (use_data & use_grad) {
      check_cudart(cudaMalloc(&trqx, rs * n));
      check_cudart(cudaMalloc(&trqy, rs * n));
      check_cudart(cudaMalloc(&trqz, rs * n));
    } else {
      trqx = nullptr;
      trqy = nullptr;
      trqz = nullptr;
    }

    check_cudart(cudaMalloc(&vir_trq, rs * 9));
  }

  if (rc & rc_copyin) {
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

  pme_data(rc);
}

void chkpole();
void rotpole();
void elec_init(int vers) {
  if (!use_elec())
    return;

  // zero torque

  if (vers & use_grad) {
    zero_array(trqx, n);
    zero_array(trqy, n);
    zero_array(trqz, n);
  }

  // zero torque-related virial

  if (vers & use_virial) {
    zero_array(vir_trq, 9);
  }

  chkpole();
  rotpole();

  if (use_ewald()) {
    pme_init(vers);
  }
}
}
TINKER_NAMESPACE_END
