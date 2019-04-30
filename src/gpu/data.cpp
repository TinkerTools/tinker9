#include "gpu/data.h"
#include <cuda_runtime.h>
#include <ext/tinker/tinker.mod.h>

TINKER_NAMESPACE_BEGIN
namespace gpu {
void copyin_data(real* dst, const double* src, int nelem) {
  size_t size = sizeof(real) * nelem;
  if (sizeof(real) == sizeof(double)) {
    cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice);
  } else if (sizeof(real) == sizeof(float)) {
    std::vector<real> buf(nelem);
    for (int i = 0; i < nelem; ++i) {
      buf[i] = src[i];
    }
    cudaMemcpy(dst, buf.data(), size, cudaMemcpyHostToDevice);
  } else {
    assert(false);
  }
}

void copyin_data3(real* d1, real* d2, real* d3, const double* src3, int nelem) {
  size_t size = sizeof(real) * nelem;
  std::vector<real> buf(nelem);
  for (int i = 0; i < nelem; ++i) {
    buf[i] = src3[3 * i];
  }
  cudaMemcpy(d1, buf.data(), size, cudaMemcpyHostToDevice);
  for (int i = 0; i < nelem; ++i) {
    buf[i] = src3[3 * i + 1];
  }
  cudaMemcpy(d2, buf.data(), size, cudaMemcpyHostToDevice);
  for (int i = 0; i < nelem; ++i) {
    buf[i] = src3[3 * i + 2];
  }
  cudaMemcpy(d3, buf.data(), size, cudaMemcpyHostToDevice);
}

int use_data = 0;

real *x, *y, *z;
void xyz_data(int op) {
  if ((use_data & use_xyz) == 0)
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
    copyin_data(x, atoms::x, n);
    copyin_data(y, atoms::y, n);
    copyin_data(z, atoms::z, n);
  }
}

real *vx, *vy, *vz;
void vel_data(int op) {
  if ((use_data & use_vel) == 0)
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
    copyin_data3(vx, vy, vz, moldyn::v, n);
  }
}

real *ax, *ay, *az;
void accel_data(int op) {
  if ((use_data & use_accel) == 0)
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
    copyin_data3(ax, ay, az, moldyn::a, n);
  }
}

real* mass;
void mass_data(int op) {
  if ((use_data & use_mass) == 0)
    return;

  if (op == op_destroy) {
    cudaFree(mass);
  }

  if (op == op_create) {
    int n = atoms::n;
    size_t size = sizeof(real) * n;
    cudaMalloc(&mass, size);
    copyin_data(mass, atomid::mass, n);
  }
}

real energy;

real *gx, *gy, *gz;
void grad_data(int op) {
  if ((use_data & use_grad) == 0)
    return;

  if (op == op_destroy) {
    cudaFree(gx);
    cudaFree(gy);
    cudaFree(gz);
  }

  if (op == op_create) {
    int n = atoms::n;
    size_t size = sizeof(real) * n;
    cudaMalloc(&gx, size);
    cudaMalloc(&gy, size);
    cudaMalloc(&gz, size);
    copyin_data3(gx, gy, gz, deriv::desum, n);
  }
}

real vir[3][3];
}

void gpu_data_create() {
  using namespace gpu;
  xyz_data(op_create);
  vel_data(op_create);
  accel_data(op_create);
  mass_data(op_create);
  grad_data(op_create);
}

void gpu_data_destroy() {
  using namespace gpu;
  xyz_data(op_destroy);
  vel_data(op_destroy);
  accel_data(op_destroy);
  mass_data(op_destroy);
  grad_data(op_destroy);
}
TINKER_NAMESPACE_END
