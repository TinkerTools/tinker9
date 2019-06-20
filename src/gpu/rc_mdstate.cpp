#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "rc_cudart.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int use_data = 0;
int n = 0;

real *x, *y, *z;
real *vx, *vy, *vz;
real *ax, *ay, *az;
real* mass;

void n_data(int op) {
  if (op & op_dealloc)
    n = 0;

  // if (op & op_alloc) {}

  if (op & op_copyin)
    n = atoms::n;
}

void xyz_data(int op) {
  if ((use_xyz & use_data) == 0)
    return;

  if (op & op_dealloc) {
    check_cudart(cudaFree(x));
    check_cudart(cudaFree(y));
    check_cudart(cudaFree(z));
  }

  if (op & op_alloc) {
    size_t size = sizeof(real) * n;
    check_cudart(cudaMalloc(&x, size));
    check_cudart(cudaMalloc(&y, size));
    check_cudart(cudaMalloc(&z, size));
  }

  if (op & op_copyin) {
    copyin_data(x, atoms::x, n);
    copyin_data(y, atoms::y, n);
    copyin_data(z, atoms::z, n);
  }
}

void vel_data(int op) {
  if ((use_vel & use_data) == 0)
    return;

  if (op & op_dealloc) {
    check_cudart(cudaFree(vx));
    check_cudart(cudaFree(vy));
    check_cudart(cudaFree(vz));
  }

  if (op & op_alloc) {
    size_t size = sizeof(real) * n;
    check_cudart(cudaMalloc(&vx, size));
    check_cudart(cudaMalloc(&vy, size));
    check_cudart(cudaMalloc(&vz, size));
  }

  if (op & op_copyin) {
    copyin_data2(0, 3, vx, moldyn::v, n);
    copyin_data2(1, 3, vy, moldyn::v, n);
    copyin_data2(2, 3, vz, moldyn::v, n);
  }
}

void accel_data(int op) {
  if ((use_accel & use_data) == 0)
    return;

  if (op & op_dealloc) {
    check_cudart(cudaFree(ax));
    check_cudart(cudaFree(ay));
    check_cudart(cudaFree(az));
  }

  if (op & op_alloc) {
    size_t size = sizeof(real) * n;
    check_cudart(cudaMalloc(&ax, size));
    check_cudart(cudaMalloc(&ay, size));
    check_cudart(cudaMalloc(&az, size));
  }

  if (op & op_copyin) {
    copyin_data2(0, 3, ax, moldyn::a, n);
    copyin_data2(1, 3, ay, moldyn::a, n);
    copyin_data2(2, 3, az, moldyn::a, n);
  }
}

void mass_data(int op) {
  if ((use_mass & use_data) == 0)
    return;

  if (op & op_dealloc) {
    check_cudart(cudaFree(mass));
  }

  if (op & op_alloc) {
    size_t size = sizeof(real) * n;
    check_cudart(cudaMalloc(&mass, size));
  }

  if (op & op_copyin) {
    copyin_data(mass, atomid::mass, n);
  }
}

real* esum;
real* vir;
real *gx, *gy, *gz;

void zero_egv() {
  int flag_e = gpu::use_data & gpu::use_energy;
  int flag_g = gpu::use_data & gpu::use_grad;
  int flag_v = gpu::use_data & gpu::use_virial;
  int n = gpu::n;

  if (flag_e) {
    gpu::zero_data(gpu::esum, 1);
  }

  if (flag_g) {
    gpu::zero_data(gpu::gx, n);
    gpu::zero_data(gpu::gy, n);
    gpu::zero_data(gpu::gz, n);
  }

  if (flag_v) {
    gpu::zero_data(gpu::vir, 9);
  }
}

void egv_data(int op) {
  if (use_mass & (use_energy + use_grad + use_virial) == 0)
    return;

  if (op & op_dealloc) {
    if (use_energy & use_data) {
      check_cudart(cudaFree(esum));
    }

    if (use_grad & use_data) {
      check_cudart(cudaFree(gx));
      check_cudart(cudaFree(gy));
      check_cudart(cudaFree(gz));
    }

    if (use_virial & use_data) {
      check_cudart(cudaFree(vir));
    }
  }

  if (op & op_alloc) {
    const size_t rs = sizeof(real);
    size_t size = 0;

    if (use_energy & use_data) {
      size = rs;
      check_cudart(cudaMalloc(&esum, size));
    }

    if (use_grad & use_data) {
      size = rs * n;
      check_cudart(cudaMalloc(&gx, size));
      check_cudart(cudaMalloc(&gy, size));
      check_cudart(cudaMalloc(&gz, size));
    }

    if (use_virial & use_data) {
      size = rs * 9;
      check_cudart(cudaMalloc(&vir, size));
    }
  }

  if (op & op_copyin) {
    if (use_energy & use_data) {
      copyin_data(esum, &energi::esum, 1);
    }

    if (use_grad & use_data) {
      copyin_data2(0, 3, gx, deriv::desum, n);
      copyin_data2(1, 3, gy, deriv::desum, n);
      copyin_data2(2, 3, gz, deriv::desum, n);
    }

    if (use_virial & use_data) {
      copyin_data(vir, &virial::vir[0][0], 9);
    }
  }
}
}
TINKER_NAMESPACE_END

extern "C" {}
