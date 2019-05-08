#include "gpu/data.h"
#include "gpu/nblist.h"
#include "util/error.cudart.h"
#include <cuda_runtime.h>

TINKER_NAMESPACE_BEGIN
namespace gpu {
void copyin_data_1(int* dst, const int* src, int nelem) {
  check_cudart(
      cudaMemcpy(dst, src, sizeof(int) * nelem, cudaMemcpyHostToDevice));
}

void copyin_data_1(real* dst, const double* src, int nelem) {
  const size_t rs = sizeof(real);
  size_t size = rs * nelem;
  if (rs == sizeof(double)) {
    check_cudart(cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice));
  } else if (rs == sizeof(float)) {
    std::vector<real> buf(nelem);
    for (int i = 0; i < nelem; ++i) {
      buf[i] = src[i];
    }
    check_cudart(cudaMemcpy(dst, buf.data(), size, cudaMemcpyHostToDevice));
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
  check_cudart(cudaMemcpy(dst, buf.data(), size, cudaMemcpyHostToDevice));
}

void n_data();
void xyz_data(int op);
void vel_data(int op);
void accel_data(int op);
void mass_data(int op);
void energy_data(int op);

void potential_data(int op);

void box_data(int op);
void couple_data(int op);
void nblist_data(int op);
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_data_create() {
  m_tinker_using_namespace;
  using namespace gpu;
  const int op = op_create;

  n_data();

  xyz_data(op);
  vel_data(op);
  accel_data(op);
  mass_data(op);
  energy_data(op);

  potential_data(op);

  box_data(op);
  couple_data(op);
  nblist_data(op);

  // build neighbor lists
  tinker_gpu_vlist_build();
}

void tinker_gpu_data_destroy() {
  m_tinker_using_namespace;
  using namespace gpu;
  const int op = op_destroy;

  xyz_data(op);
  vel_data(op);
  accel_data(op);
  mass_data(op);
  energy_data(op);

  potential_data(op);

  // Neighbor lists must be initialized after potential initialization.
  // xred, yred, and zred need to be initialized in vdw routines and will be
  // used in nblist setups.
  box_data(op);
  couple_data(op);
  nblist_data(op);
}
}
