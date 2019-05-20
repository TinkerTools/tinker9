#ifndef TINKER_GPU_ACC_FMAT_H_
#define TINKER_GPU_ACC_FMAT_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
/**
 * @brief
 * A Fortran matrix f(5,3) is equivalent to a 2D C array c[3][5], with NFRow = 5
 * and NFCol = 3.
 */
template <class T, int NFRow, int NFCol>
class FortranMatrixView {
private:
  T (&data_)[NFCol][NFRow];

public:
  #pragma acc routine seq
  FortranMatrixView(T (&_data)[NFCol][NFRow]) : data_(_data) {}

  #pragma acc routine seq
  const T& operator()(int ifrow1, int ifcol1) const {
    return data_[ifcol1 - 1][ifrow1 - 1];
  }

  #pragma acc routine seq
  T& operator()(int ifrow1, int ifcol1) {
    return data_[ifcol1 - 1][ifrow1 - 1];
  }
};

template <int NFRow, int NFCol>
using freal = FortranMatrixView<real, NFRow, NFCol>;
using freal33 = FortranMatrixView<real, 3, 3>;
}
TINKER_NAMESPACE_END

#endif
