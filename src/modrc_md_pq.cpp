#include "mod_box.h"
#include "mod_md.h"
#include "util_array.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN

//======================================================================
// number of atoms

void n_data(rc_op op) {
  if (op & rc_dealloc) {
    trajn = -1;
    n = 0;
  }

  if (op & rc_alloc) {
    n = atoms::n;

    if (calc::traj & use_data) {
      // trajn must have been initialized by this point
      assert(trajn >= 0);
    }
  }
}

//======================================================================
// x y z coordinates

void xyz_data(rc_op op) {
  if ((calc::xyz & use_data) == 0)
    return;

  if (op & rc_dealloc) {
    if (calc::traj & use_data) {
      check_rt(cudaFree(trajx));
      check_rt(cudaFree(trajy));
      check_rt(cudaFree(trajz));
      x = nullptr;
      y = nullptr;
      z = nullptr;
    } else {
      trajx = nullptr;
      trajy = nullptr;
      trajz = nullptr;
      check_rt(cudaFree(x));
      check_rt(cudaFree(y));
      check_rt(cudaFree(z));
    }
  }

  if (op & rc_alloc) {
    size_t size = sizeof(real) * n;
    if (calc::traj & use_data) {
      size *= trajn;
      check_rt(cudaMalloc(&trajx, size));
      check_rt(cudaMalloc(&trajy, size));
      check_rt(cudaMalloc(&trajz, size));
      x = trajx;
      y = trajy;
      z = trajz;
    } else {
      check_rt(cudaMalloc(&x, size));
      check_rt(cudaMalloc(&y, size));
      check_rt(cudaMalloc(&z, size));
    }
  }

  if (op & rc_init) {
    copyin_array(x, atoms::x, n);
    copyin_array(y, atoms::y, n);
    copyin_array(z, atoms::z, n);
  }
}

//======================================================================
// velocitiies

void vel_data(rc_op op) {
  if ((calc::vel & use_data) == 0)
    return;

  if (op & rc_dealloc) {
    check_rt(cudaFree(vx));
    check_rt(cudaFree(vy));
    check_rt(cudaFree(vz));
  }

  if (op & rc_alloc) {
    size_t size = sizeof(real) * n;
    check_rt(cudaMalloc(&vx, size));
    check_rt(cudaMalloc(&vy, size));
    check_rt(cudaMalloc(&vz, size));
  }

  if (op & rc_init) {
    copyin_array2(0, 3, vx, moldyn::v, n);
    copyin_array2(1, 3, vy, moldyn::v, n);
    copyin_array2(2, 3, vz, moldyn::v, n);
  }
}

//======================================================================
// atomic mass

void mass_data(rc_op op) {
  if ((calc::mass & use_data) == 0)
    return;

  if (op & rc_dealloc) {
    check_rt(cudaFree(mass));
    check_rt(cudaFree(massinv));
  }

  if (op & rc_alloc) {
    size_t size = sizeof(real) * n;
    check_rt(cudaMalloc(&mass, size));
    check_rt(cudaMalloc(&massinv, size));
  }

  if (op & rc_init) {
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
