#include "gpu/decl_basic.h"
#include "gpu/rc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void zero_array(real* dst, int nelem) {
  size_t size = sizeof(real) * nelem;
  check_cudart(cudaMemset(dst, 0, size));
}

void copyin_array(int* dst, const int* src, int nelem) {
  check_cudart(
      cudaMemcpy(dst, src, sizeof(int) * nelem, cudaMemcpyHostToDevice));
}

void copyin_array(real* dst, const double* src, int nelem) {
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

void copyin_array2(int idx0, int ndim, real* dst, const double* src,
                   int nelem) {
  size_t size = sizeof(real) * nelem;
  std::vector<real> buf(nelem);
  for (int i = 0; i < nelem; ++i) {
    buf[i] = src[ndim * i + idx0];
  }
  check_cudart(cudaMemcpy(dst, buf.data(), size, cudaMemcpyHostToDevice));
}

void copyout_array(int* dst, const int* src, int nelem) {
  check_cudart(
      cudaMemcpy(dst, src, sizeof(int) * nelem, cudaMemcpyDeviceToHost));
}

void copyout_array(double* dst, const real* src, int nelem) {
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

void copyout_array2(int idx0, int ndim, double* dst, const real* src,
                    int nelem) {
  size_t size = sizeof(real) * nelem;
  std::vector<real> buf(nelem);
  check_cudart(cudaMemcpy(buf.data(), src, size, cudaMemcpyDeviceToHost));
  for (int i = 0; i < nelem; ++i) {
    dst[ndim * i + idx0] = buf[i];
  }
}

void copyout_array3(double (*dst)[3], const real (*src)[3], int natom) {
  copyout_array(&dst[0][0], &src[0][0], 3 * natom);
}

void copyout_array3(std::vector<std::array<double, 3>>& dst,
                    const real (*src)[3], int natom) {
  dst.resize(natom);
  copyout_array(&dst[0][0], &src[0][0], 3 * natom);
}

void copy_array(int* dst, const int* src, int nelem) {
  size_t size = sizeof(int) * nelem;
  check_cudart(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToDevice));
}

void copy_array(real* dst, const real* src, int nelem) {
  size_t size = sizeof(real) * nelem;
  check_cudart(cudaMemcpy(dst, src, size, cudaMemcpyDeviceToDevice));
}

}
TINKER_NAMESPACE_END
