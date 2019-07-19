#include "gpu/decl_mdstate.h"
#include "gpu/rc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int use_data;

//======================================================================
// number of atoms

int trajn = -1;
int n;

void n_data(rc_t rc) {
  if (rc & rc_dealloc) {
    trajn = -1;
    n = 0;
  }

  if (rc & rc_alloc) {
    n = atoms::n;

    if (use_traj & use_data) {
      // trajn must have been initialized by this point
      assert(trajn >= 0);
    }
  }
}

//======================================================================
// x y z coordinates

real *trajx, *trajy, *trajz;
real *x, *y, *z;

void xyz_data(rc_t rc) {
  if ((use_xyz & use_data) == 0)
    return;

  if (rc & rc_dealloc) {
    if (use_traj & use_data) {
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
    if (use_traj & use_data) {
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

real *vx, *vy, *vz;

void vel_data(rc_t rc) {
  if ((use_vel & use_data) == 0)
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

real* mass;
real* massinv;

void mass_data(rc_t rc) {
  if ((use_mass & use_data) == 0)
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

extern void potential_data_(rc_t);
extern void md_data(rc_t);
void mdstate_data(rc_t rc) {
  n_data(rc);

  xyz_data(rc);
  vel_data(rc);
  // accel_data(rc);
  mass_data(rc);

  potential_data_(rc);

  // Neighbor lists must be initialized after potential initialization.
  // xred, yred, and zred need to be initialized in vdw routines and will be
  // used in nblist setups.
  box_data(rc);
  couple_data(rc);
  nblist_data(rc);

  random_data(rc);

  md_data(rc);
}

void goto_frame0(int idx0) {
  assert(use_traj & use_data);
  x = trajx + n * idx0;
  y = trajy + n * idx0;
  z = trajz + n * idx0;
  box = trajbox + idx0;
}

void goto_frame1(int idx1) { goto_frame0(idx1 - 1); }
}
TINKER_NAMESPACE_END
