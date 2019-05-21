#ifndef TINKER_GPU_ACC_FMAT_H_
#define TINKER_GPU_ACC_FMAT_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
/**
 * @brief
 * A Fortran matrix f(3,50) is equivalent to a 2D C array c[50][3],
 * with NFRow = 3 and NFCol = 50.
 */
template <class T, int NFRow>
class FortranMatrixView {
private:
  T (*data_)[NFRow];

public:
  #pragma acc routine seq
  FortranMatrixView(const T (*_data)[NFRow])
      : data_(const_cast<T (*)[NFRow]>(_data)) {}
  #pragma acc routine seq
  FortranMatrixView(T (*_data)[NFRow]) : data_(_data) {}

  #pragma acc routine seq
  const T& operator()(int ifrow1, int ifcol1) const {
    return data_[ifcol1 - 1][ifrow1 - 1];
  }

  #pragma acc routine seq
  T& operator()(int ifrow1, int ifcol1) {
    return data_[ifcol1 - 1][ifrow1 - 1];
  }
};

template <int NFRow>
using freal_mat = FortranMatrixView<real, NFRow>;
using freal_mat3 = FortranMatrixView<real, 3>;
}
TINKER_NAMESPACE_END

#endif
