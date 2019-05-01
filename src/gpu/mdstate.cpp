#include "gpu/mdstate.h"
#include "util/cxx.h"
#include <cuda_runtime.h>
#include <ext/tinker/tinker.mod.h>

TINKER_NAMESPACE_BEGIN
namespace gpu {
void copyin_data_1(int* dst, const int* src, int nelem) {
  cudaMemcpy(dst, src, sizeof(int) * nelem, cudaMemcpyHostToDevice);
}

void copyin_data_1(real* dst, const double* src, int nelem) {
  const size_t rs = sizeof(real);
  size_t size = rs * nelem;
  if (rs == sizeof(double)) {
    cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice);
  } else if (rs == sizeof(float)) {
    std::vector<real> buf(nelem);
    for (int i = 0; i < nelem; ++i) {
      buf[i] = src[i];
    }
    cudaMemcpy(dst, buf.data(), size, cudaMemcpyHostToDevice);
  } else {
    assert(false);
  }
}

void copyin_data_n(int idx0, int ndim, real* dst, const double* src,
                   int nelem) {
  size_t size = sizeof(real) * nelem;
  std::vector<real> buf(nelem);
  for (int i = 0; i < nelem; ++i) {
    buf[i] = src[ndim * i + idx0];
  }
  cudaMemcpy(dst, buf.data(), size, cudaMemcpyHostToDevice);
}

int use_data = 0;

real* esum;
real* vir;
real* mass;

real *x, *y, *z;
real *vx, *vy, *vz;
real *ax, *ay, *az;
real *gx, *gy, *gz;

void xyz_data(int op) {
  if ((use_xyz & use_data) == 0)
    return;

  if (op == op_destroy) {
    cudaFree(x);
    cudaFree(y);
    cudaFree(z);
  }

  if (op == op_create) {
    int n = atoms::n;
    size_t size = sizeof(real) * n;
    cudaMalloc(&x, size);
    cudaMalloc(&y, size);
    cudaMalloc(&z, size);
    copyin_data_1(x, atoms::x, n);
    copyin_data_1(y, atoms::y, n);
    copyin_data_1(z, atoms::z, n);
  }
}

void vel_data(int op) {
  if ((use_vel & use_data) == 0)
    return;

  if (op == op_destroy) {
    cudaFree(vx);
    cudaFree(vy);
    cudaFree(vz);
  }

  if (op == op_create) {
    int n = atoms::n;
    size_t size = sizeof(real) * n;
    cudaMalloc(&vx, size);
    cudaMalloc(&vy, size);
    cudaMalloc(&vz, size);
    copyin_data_n(0, 3, vx, moldyn::v, n);
    copyin_data_n(1, 3, vy, moldyn::v, n);
    copyin_data_n(2, 3, vz, moldyn::v, n);
  }
}

void accel_data(int op) {
  if ((use_accel & use_data) == 0)
    return;

  if (op == op_destroy) {
    cudaFree(ax);
    cudaFree(ay);
    cudaFree(az);
  }

  if (op == op_create) {
    int n = atoms::n;
    size_t size = sizeof(real) * n;
    cudaMalloc(&ax, size);
    cudaMalloc(&ay, size);
    cudaMalloc(&az, size);
    copyin_data_n(0, 3, ax, moldyn::a, n);
    copyin_data_n(1, 3, ay, moldyn::a, n);
    copyin_data_n(2, 3, az, moldyn::a, n);
  }
}

void mass_data(int op) {
  if ((use_mass & use_data) == 0)
    return;

  if (op == op_destroy) {
    cudaFree(mass);
  }

  if (op == op_create) {
    int n = atoms::n;
    size_t size = sizeof(real) * n;
    cudaMalloc(&mass, size);
    copyin_data_1(mass, atomid::mass, n);
  }
}

void energy_data(int op) {
  if (use_mass & (use_energy + use_grad + use_virial) == 0)
    return;

  if (op == op_destroy) {
    if (use_energy & use_data) {
      cudaFree(esum);
    }

    if (use_grad & use_data) {
      cudaFree(gx);
      cudaFree(gy);
      cudaFree(gz);
    }

    if (use_virial & use_data) {
      cudaFree(vir);
    }
  }

  if (op == op_create) {
    const size_t rs = sizeof(real);
    size_t size = 0;
    if (use_energy & use_data) {
      size = rs;
      cudaMalloc(&esum, size);
      copyin_data_1(esum, &energi::esum, 1);
    }

    if (use_grad & use_data) {
      int n = atoms::n;
      size = rs * n;
      cudaMalloc(&gx, size);
      cudaMalloc(&gy, size);
      cudaMalloc(&gz, size);
      copyin_data_n(0, 3, gx, deriv::desum, n);
      copyin_data_n(1, 3, gy, deriv::desum, n);
      copyin_data_n(2, 3, gz, deriv::desum, n);
    }

    if (use_virial & use_data) {
      size = rs * 9;
      cudaMalloc(&vir, size);
      copyin_data_1(vir, &virial::vir[0][0], 9);
    }
  }
}

box_st* box;

void box_data(int op) {
  if (op == op_destroy) {
    cudaFree(box);
  }

  if (op == op_create) {
    int shape = box_null;
    if (boxes::orthogonal)
      shape = box_ortho;
    else if (boxes::monoclinic)
      shape = box_mono;
    else if (boxes::triclinic)
      shape = box_tri;
    else if (boxes::octahedron)
      shape = box_oct;

    size_t size = sizeof(box_st);
    cudaMalloc(&box, size);

    copyin_data_1(&box->xbox, &boxes::xbox, 1);
    copyin_data_1(&box->ybox, &boxes::ybox, 1);
    copyin_data_1(&box->zbox, &boxes::zbox, 1);
    copyin_data_1(&box->alpha, &boxes::alpha, 1);
    copyin_data_1(&box->beta, &boxes::beta, 1);
    copyin_data_1(&box->gamma, &boxes::gamma, 1);
    copyin_data_1(&box->lvec[0][0], &boxes::lvec[0][0], 9);
    copyin_data_1(&box->recip[0][0], &boxes::recip[0][0], 9);
    copyin_data_1(&box->shape, &shape, 1);
  }
}
}
TINKER_NAMESPACE_END
