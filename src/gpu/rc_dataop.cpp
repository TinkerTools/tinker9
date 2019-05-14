#include "gpu/decl_dataop.h"
#include "rc_cudart.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void copyin_data(int* dst, const int* src, int nelem) {
  check_cudart(
      cudaMemcpy(dst, src, sizeof(int) * nelem, cudaMemcpyHostToDevice));
}

void copyin_data(real* dst, const double* src, int nelem) {
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

void copyin_data2(int idx0, int ndim, real* dst, const double* src, int nelem) {
  size_t size = sizeof(real) * nelem;
  std::vector<real> buf(nelem);
  for (int i = 0; i < nelem; ++i) {
    buf[i] = src[ndim * i + idx0];
  }
  check_cudart(cudaMemcpy(dst, buf.data(), size, cudaMemcpyHostToDevice));
}

void copyout_data(int* dst, const int* src, int nelem) {
  check_cudart(
      cudaMemcpy(dst, src, sizeof(int) * nelem, cudaMemcpyDeviceToHost));
}

void copyout_data(double* dst, const real* src, int nelem) {
  const size_t rs = sizeof(real);
  size_t size = rs * nelem;
  if (rs == sizeof(double)) {
    check_cudart(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToHost));
  } else if (rs == sizeof(float)) {
    std::vector<real> buf(nelem);
    check_cudart(cudaMemcpy(buf.data(), src, size, cudaMemcpyDeviceToHost));
    for (int i = 0; i < nelem; ++i) {
      dst[i] = buf[i];
    }
  } else {
    assert(false);
  }
}

void copyout_data2(int idx0, int ndim, double* dst, const real* src,
                   int nelem) {
  size_t size = sizeof(real) * nelem;
  std::vector<real> buf(nelem);
  check_cudart(cudaMemcpy(buf.data(), src, size, cudaMemcpyDeviceToHost));
  for (int i = 0; i < nelem; ++i) {
    dst[ndim * i + idx0] = buf[i];
  }
}

void zero_data(real* dst, int nelem) {
  size_t size = sizeof(real) * nelem;
  check_cudart(cudaMemset(dst, 0, size));
}
}
TINKER_NAMESPACE_END

#include "gpu/decl_box.h"
#include "gpu/decl_couple.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_nblist.h"
#include "gpu/e_potential.h"

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
  egv_data(op);

  potential_data(op);

  // Neighbor lists must be initialized after potential initialization.
  // xred, yred, and zred need to be initialized in vdw routines and will be
  // used in nblist setups.
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
  egv_data(op);

  potential_data(op);

  box_data(op);
  couple_data(op);
  nblist_data(op);
}
}
