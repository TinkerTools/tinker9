#include "gpu/decl_mdstate.h"
#include "gpu/rc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int use_data;

//======================================================================
// number of atoms

int n;

void n_data(rc_t rc) {
  if (rc & rc_dealloc) {
    n = 0;
  }

  if (rc & rc_alloc) {
    n = atoms::n;
  }

  if (rc & rc_copyin) {
    n = atoms::n;
  }
}

//======================================================================
/// x y z coordinates
real *x, *y, *z;

void xyz_data(rc_t rc) {
  if ((use_xyz & use_data) == 0)
    return;

  if (rc & rc_dealloc) {
    check_cudart(cudaFree(x));
    check_cudart(cudaFree(y));
    check_cudart(cudaFree(z));
  }

  if (rc & rc_alloc) {
    size_t size = sizeof(real) * n;
    check_cudart(cudaMalloc(&x, size));
    check_cudart(cudaMalloc(&y, size));
    check_cudart(cudaMalloc(&z, size));
  }

  if (rc & rc_copyin) {
    copyin_array(x, atoms::x, n);
    copyin_array(y, atoms::y, n);
    copyin_array(z, atoms::z, n);
  }

  if (rc & rc_copyout) {
    copyout_array(atoms::x, x, n);
    copyout_array(atoms::y, y, n);
    copyout_array(atoms::z, z, n);
  }
}

//======================================================================
/// velocitiies
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

  if (rc & rc_copyout) {
    copyout_array2(0, 3, moldyn::v, vx, n);
    copyout_array2(1, 3, moldyn::v, vy, n);
    copyout_array2(2, 3, moldyn::v, vz, n);
  }
}

//======================================================================
/// accelerations
real *ax, *ay, *az;

void accel_data(rc_t rc) {
  if ((use_accel & use_data) == 0)
    return;

  if (rc & rc_dealloc) {
    check_cudart(cudaFree(ax));
    check_cudart(cudaFree(ay));
    check_cudart(cudaFree(az));
  }

  if (rc & rc_alloc) {
    size_t size = sizeof(real) * n;
    check_cudart(cudaMalloc(&ax, size));
    check_cudart(cudaMalloc(&ay, size));
    check_cudart(cudaMalloc(&az, size));
  }

  if (rc & rc_copyin) {
    copyin_array2(0, 3, ax, moldyn::a, n);
    copyin_array2(1, 3, ay, moldyn::a, n);
    copyin_array2(2, 3, az, moldyn::a, n);
  }

  if (rc & rc_copyout) {
    copyout_array2(0, 3, moldyn::a, ax, n);
    copyout_array2(1, 3, moldyn::a, ay, n);
    copyout_array2(2, 3, moldyn::a, az, n);
  }
}

//======================================================================
/// atomic mass
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

  // if (rc & rc_copyout)
  //   copyout_array(atomid::mass, mass, n);
}

extern void potential_data__(rc_t);
void mdstate_data(rc_t rc) {
  n_data(rc);

  xyz_data(rc);
  vel_data(rc);
  accel_data(rc);
  mass_data(rc);

  potential_data__(rc);

  // Neighbor lists must be initialized after potential initialization.
  // xred, yred, and zred need to be initialized in vdw routines and will be
  // used in nblist setups.
  box_data(rc);
  couple_data(rc);
  nblist_data(rc);
}
}
TINKER_NAMESPACE_END
