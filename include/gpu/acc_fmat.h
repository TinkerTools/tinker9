#ifndef TINKER_GPU_ACC_FMAT_H_
#define TINKER_GPU_ACC_FMAT_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <class T, int... Ints>
class FortranMatrixView;

template <class T>
class FortranMatrixView<T> {
private:
  T* data_;

public:
  #pragma acc routine seq
  FortranMatrixView(const T* _data) : data_(const_cast<T*>(_data)) {}

  #pragma acc routine seq
  FortranMatrixView(T* _data) : data_(_data) {}

  #pragma acc routine seq
  const T& operator()(int if1) const { return data_[if1 - 1]; }

  #pragma acc routine seq
  T& operator()(int if1) { return data_[if1 - 1]; }
};

using farray_real = FortranMatrixView<real>;

/**
 * @brief
 * A Fortran matrix f(3,50) is equivalent to a 2D C array c[50][3],
 * with NFRow = 3 and NFCol = 50.
 */
template <class T, int NFRow>
class FortranMatrixView<T, NFRow> {
private:
  T* data_;

public:
  #pragma acc routine seq
  FortranMatrixView(const T* _data) : data_(const_cast<T*>(_data)) {}

  #pragma acc routine seq
  FortranMatrixView(T* _data) : data_(_data) {}

  #pragma acc routine seq
  FortranMatrixView(const T (*_data)[NFRow])
      : data_(const_cast<T*>(&_data[0][0])) {}

  #pragma acc routine seq
  FortranMatrixView(T (*_data)[NFRow]) : data_(&_data[0][0]) {}

  #pragma acc routine seq
  const T& operator()(int ifrow1, int ifcol1) const {
    return data_[(ifcol1 - 1) * NFRow + ifrow1 - 1];
  }

  #pragma acc routine seq
  T& operator()(int ifrow1, int ifcol1) {
    return data_[(ifcol1 - 1) * NFRow + ifrow1 - 1];
  }
};

template <int NFRow>
using fmat_real = FortranMatrixView<real, NFRow>;
using fmat_real3 = FortranMatrixView<real, 3>;

template <class T>
class FortranDynamicMatrixView {
private:
  T* data_;
  int nfrow_;

public:
  #pragma acc routine seq
  FortranDynamicMatrixView(const T* _data, int _dim)
      : data_(const_cast<T*>(_data)), nfrow_(_dim) {}

  #pragma acc routine seq
  FortranDynamicMatrixView(T* _data, int _dim) : data_(_data), nfrow_(_dim) {}

  #pragma acc routine seq
  FortranDynamicMatrixView(const T* _data, int _dim, int /* dummy_dim */)
      : data_(const_cast<T*>(_data)), nfrow_(_dim) {}

  #pragma acc routine seq
  FortranDynamicMatrixView(T* _data, int _dim, int /* dummy_dim */)
      : data_(_data), nfrow_(_dim) {}

  #pragma acc routine seq
  const T& operator()(int ifrow1, int ifcol1) const {
    return data_[(ifcol1 - 1) * nfrow_ + ifrow1 - 1];
  }

  #pragma acc routine seq
  T& operator()(int ifrow1, int ifcol1) {
    return data_[(ifcol1 - 1) * nfrow_ + ifrow1 - 1];
  }
};

template <class T>
using allocatable = FortranDynamicMatrixView<T>;
using real_allocatable = allocatable<real>;
}
TINKER_NAMESPACE_END

#endif
