#include "mod_box.h"
#include "mod_md.h"
#include "util_array.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN

//======================================================================
// number of atoms

void n_data(rc_t rc) {
  if (rc & rc_dealloc) {
    trajn = -1;
    n = 0;
  }

  if (rc & rc_alloc) {
    n = atoms::n;

    if (calc::traj & use_data) {
      // trajn must have been initialized by this point
      assert(trajn >= 0);
    }
  }
}

//======================================================================
// x y z coordinates

void xyz_data(rc_t rc) {
  if ((calc::xyz & use_data) == 0)
    return;

  if (rc & rc_dealloc) {
    if (calc::traj & use_data) {
      check_cudart(cudaFree(trajx));
      check_cudart(cudaFree(trajy));
      check_cudart(cudaFree(trajz));
      x = nullptr;
      y = nullptr;
      z = nullptr;
    } else {
      trajx = nullptr;
      trajy = nullptr;
      trajz = nullptr;
      check_cudart(cudaFree(x));
      check_cudart(cudaFree(y));
      check_cudart(cudaFree(z));
    }
  }

  if (rc & rc_alloc) {
    size_t size = sizeof(real) * n;
    if (calc::traj & use_data) {
      size *= trajn;
      check_cudart(cudaMalloc(&trajx, size));
      check_cudart(cudaMalloc(&trajy, size));
      check_cudart(cudaMalloc(&trajz, size));
      x = trajx;
      y = trajy;
      z = trajz;
    } else {
      check_cudart(cudaMalloc(&x, size));
      check_cudart(cudaMalloc(&y, size));
      check_cudart(cudaMalloc(&z, size));
    }
  }

  if (rc & rc_copyin) {
    copyin_array(x, atoms::x, n);
    copyin_array(y, atoms::y, n);
    copyin_array(z, atoms::z, n);
  }
}

//======================================================================
// velocitiies

void vel_data(rc_t rc) {
  if ((calc::vel & use_data) == 0)
    return;

  if (rc & rc_dealloc) {
    check_cudart(cudaFree(vx));
    check_cudart(cudaFree(vy));
    check_cudart(cudaFree(vz));
  }

  if (rc & rc_alloc) {
    size_t size = sizeof(real) * n;
    check_cudart(cudaMalloc(&vx, size));
    check_cudart(cudaMalloc(&vy, size));
    check_cudart(cudaMalloc(&vz, size));
  }

  if (rc & rc_copyin) {
    copyin_array2(0, 3, vx, moldyn::v, n);
    copyin_array2(1, 3, vy, moldyn::v, n);
    copyin_array2(2, 3, vz, moldyn::v, n);
  }
}

//======================================================================
// atomic mass

void mass_data(rc_t rc) {
  if ((calc::mass & use_data) == 0)
    return;

  if (rc & rc_dealloc) {
    check_cudart(cudaFree(mass));
    check_cudart(cudaFree(massinv));
  }

  if (rc & rc_alloc) {
    size_t size = sizeof(real) * n;
    check_cudart(cudaMalloc(&mass, size));
    check_cudart(cudaMalloc(&massinv, size));
  }

  if (rc & rc_copyin) {
    copyin_array(mass, atomid::mass, n);
    std::vector<double> mbuf(n);
    for (int i = 0; i < n; ++i)
      mbuf[i] = 1 / atomid::mass[i];
    copyin_array(massinv, mbuf.data(), n);
  }
}

void goto_frame0(int idx0) {
  assert(calc::traj & use_data);
  x = trajx + n * idx0;
  y = trajy + n * idx0;
  z = trajz + n * idx0;
  box = trajbox + idx0;
}

void goto_frame1(int idx1) { goto_frame0(idx1 - 1); }

TINKER_NAMESPACE_END
