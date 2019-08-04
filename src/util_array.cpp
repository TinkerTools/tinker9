#include "util_array.h"
#include "util_rt.h"

TINKER_NAMESPACE_BEGIN
template <class T>
void zero_array_tmpl(T* dst, int nelem) {
  size_t size = sizeof(T) * nelem;
  zero_bytes(dst, size);
}

void zero_array(int* dst, int nelem) { zero_array_tmpl(dst, nelem); }

void zero_array(float* dst, int nelem) { zero_array_tmpl(dst, nelem); }

void zero_array(double* dst, int nelem) { zero_array_tmpl(dst, nelem); }

template <class DT, class ST>
void copyin_array_tmpl(DT* dst, const ST* src, int nelem) {
  constexpr size_t ds = sizeof(DT);
  constexpr size_t ss = sizeof(ST);
  static_assert(ds <= ss, "invalid if dst = double and src = float");

  size_t size = ds * nelem;
  if_constexpr(ds == ss) { copyin_bytes(dst, src, size); }
  else if_constexpr(ds < ss) {
    std::vector<DT> buf(nelem);
    for (int i = 0; i < nelem; ++i)
      buf[i] = src[i];
    copyin_bytes(dst, buf.data(), size);
  }
}

template <class DT, class ST>
void copyout_array_tmpl(DT* dst, const ST* src, int nelem) {
  constexpr size_t ds = sizeof(DT);
  constexpr size_t ss = sizeof(ST);
  static_assert(ds >= ss, "invalid if dst = float and src = double");

  size_t size = ss * nelem;
  if_constexpr(ds == ss) { copyout_bytes(dst, src, size); }
  else if_constexpr(ds > ss) {
    std::vector<ST> buf(nelem);
    copyout_bytes(buf.data(), src, size);
    for (int i = 0; i < nelem; ++i)
      dst[i] = buf[i];
  }
}

void copyin_array(int* dst, const int* src, int nelem) {
  copyin_array_tmpl(dst, src, nelem);
}

void copyout_array(int* dst, const int* src, int nelem) {
  copyout_array_tmpl(dst, src, nelem);
}

void copyin_array(float* dst, const float* src, int nelem) {
  copyin_array_tmpl(dst, src, nelem);
}

void copyout_array(float* dst, const float* src, int nelem) {
  copyout_array_tmpl(dst, src, nelem);
}

void copyin_array(float* dst, const double* src, int nelem) {
  copyin_array_tmpl(dst, src, nelem);
}

void copyout_array(double* dst, const float* src, int nelem) {
  copyout_array_tmpl(dst, src, nelem);
}

void copyin_array(double* dst, const double* src, int nelem) {
  copyin_array_tmpl(dst, src, nelem);
}

void copyout_array(double* dst, const double* src, int nelem) {
  copyout_array_tmpl(dst, src, nelem);
}

template <class DT, class ST>
void copyin_array2_tmpl(int idx0, int ndim, DT* dst, const ST* src, int nelem) {
  std::vector<DT> buf(nelem);
  for (int i = 0; i < nelem; ++i)
    buf[i] = src[ndim * i + idx0];
  copyin_array(dst, buf.data(), nelem);
}

template <class DT, class ST>
void copyout_array2_tmpl(int idx0, int ndim, DT* dst, const ST* src,
                         int nelem) {
  std::vector<ST> buf(nelem);
  copyout_array(buf.data(), src, nelem);
  for (int i = 0; i < nelem; ++i)
    dst[ndim * i + idx0] = buf[i];
}

void copyin_array2(int idx0, int ndim, float* dst, const float* src,
                   int nelem) {
  copyin_array2_tmpl(idx0, ndim, dst, src, nelem);
}
void copyout_array2(int idx0, int ndim, float* dst, const float* src,
                    int nelem) {
  copyout_array2_tmpl(idx0, ndim, dst, src, nelem);
}

void copyin_array2(int idx0, int ndim, float* dst, const double* src,
                   int nelem) {
  copyin_array2_tmpl(idx0, ndim, dst, src, nelem);
}
void copyout_array2(int idx0, int ndim, double* dst, const float* src,
                    int nelem) {
  copyout_array2_tmpl(idx0, ndim, dst, src, nelem);
}

void copyin_array2(int idx0, int ndim, double* dst, const double* src,
                   int nelem) {
  copyin_array2_tmpl(idx0, ndim, dst, src, nelem);
}
void copyout_array2(int idx0, int ndim, double* dst, const double* src,
                    int nelem) {
  copyout_array2_tmpl(idx0, ndim, dst, src, nelem);
}

template <class DT, class ST>
void copy_array_tmpl(DT* dst, const ST* src, int nelem) {
  static_assert(std::is_same<DT, ST>::value, "");
  size_t size = sizeof(ST) * nelem;
  copy_bytes(dst, src, size);
}

void copy_array(int* dst, const int* src, int nelem) {
  copy_array_tmpl(dst, src, nelem);
}

void copy_array(float* dst, const float* src, int nelem) {
  copy_array_tmpl(dst, src, nelem);
}

void copy_array(double* dst, const double* src, int nelem) {
  copy_array_tmpl(dst, src, nelem);
}
TINKER_NAMESPACE_END
