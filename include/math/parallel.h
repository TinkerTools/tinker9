#pragma once
#include "math/parallelacc.h"
#include "math/parallelcu.h"
#include "tool/externfunc.h"
#include "tool/macro.h"

namespace tinker {
/// \ingroup math_parallel
/// \brief Sum over all of the elements of an 1D array.
///
/// \f[ Sum = \sum_i^n a_i \f]
/// \return The sum.
template <class T>
T reduceSum(const T* gpu_a, size_t nelem, int queue)
{
   // TINKER_FVOID2(cu, 1, acc, 1, T, reduceSum, const T*, size_t, int);
   return TINKER_FCALL2(cu, 1, acc, 1, reduceSum, gpu_a, nelem, queue);
}

/// \ingroup math_parallel
/// \brief Sum over all of the elements of a 2D array.
///
/// In Fortran syntax:
/// \f[ Ans(k) = \sum_i^n v(k,i), 1 \le k \le HN \f]
/// In C++ syntax:
/// \f[ Ans[k] = \sum_i^n v[i][k], 0 \le k < HN \f]
template <class HT, size_t HN, class DPTR>
void reduceSum2(HT (&h_ans)[HN], DPTR v, size_t nelem, int queue)
{
   // TINKER_FVOID2(cu, 1, acc, 1, HT (&)[HN], DPTR, size_t, int);
   TINKER_FCALL2(cu, 1, acc, 1, reduceSum2, h_ans, v, nelem, queue);
}

/// \ingroup math_parallel
/// \brief Sum over all of the elements of an 1D array. This routine will save
/// the result on the device memory in an asynchronous/non-blocking manner.
/// A valid device pointer is required for the result on the device.
///
/// \f[ Sum = \sum_i^n a_i \f]
///
/// \param dp_ans  Device pointer used to store the reduction result.
/// \param a       Device pointer to the array.
/// \param nelem   Number of elements.
/// \param queue   OpenACC queue.
template <class T>
void reduceSumOnDevice(T* dp_ans, const T* a, size_t nelem, int queue)
{
   // TINKER_FVOID2(cu, 1, acc, 1, reduceSumOnDevice, T*, const T*, size_t, int);
   TINKER_FCALL2(cu, 1, acc, 1, reduceSumOnDevice, dp_ans, a, nelem, queue);
}

/// \ingroup math_parallel
/// \brief Sum over all of the elements of a 2D array. This routine will save
/// the result on the device memory in an asynchronous/non-blocking manner.
/// A valid device array is required for the result on device.
///
/// Fortran syntax:
/// \f[ Ans(k) = \sum_i^n v(k,i), 1 \le k \le HN \f]
/// C++ syntax:
/// \f[ Ans[k] = \sum_i^n v[i][k], 0 \le k < HN \f]
///
/// \tparam HT    Type of the array element.
/// \tparam HN    Length of the result array.
/// \tparam DPTR  Type of the 2D array.
/// \param dref   Reference to the device array that stores the reduction result.
/// \param v      Device pointer to the 2D array.
/// \param nelem  Number of elements.
/// \param queue  OpenACC queue.
template <class HT, size_t HN, class DPTR>
void reduceSum2OnDevice(HT (&dref)[HN], DPTR v, size_t nelem, int queue)
{
   // TINKER_FVOID2(cu, 1, acc, 1 reduceSum2OnDevice, HT(&)[HN], DPTR, size_t, int);
   TINKER_FCALL2(cu, 1, acc, 1, reduceSum2OnDevice, dref, v, nelem, queue);
}

/// \ingroup math_parallel
/// \brief Dot product of two linear arrays.
///
/// \f[ DotProduct = \sum_i^n a_i \cdot b_i \f]
/// \return The dot product to the host thread.
template <class T>
T dotProd(const T* a, const T* b, size_t nelem, int queue)
{
   // TINKER_FVOID1(cu, 0, acc, 1, T, dotProd, const T*, const T*, size_t, int);
   return TINKER_FCALL1(cu, 0, acc, 1, dotProd, a, b, nelem, queue);
}

/// \ingroup math_parallel
/// \brief Dot product of two linear arrays.
template <class T>
void dotProd(T* ans, const T* a, const T* b, size_t nelem, int queue)
{
   // TINKER_FVOID2(cu, 1, acc, 1, dotProd, T*, const T*, const T*, size_t, int);
   TINKER_FCALL2(cu, 1, acc, 1, dotProd, ans, a, b, nelem, queue);
}

/// \ingroup math_parallel
/// \brief Multiply all of the elements in an 1D array by a scalar.
///
/// \f[ a_i = c \cdot a_i \f]
template <class T>
void scaleArray(T* dst, T scal, size_t nelem, int queue)
{
   // TINKER_FVOID2(cu, 1, acc, 1, scaleArray, T*, T, size_t, int);
   TINKER_FCALL2(cu, 1, acc, 1, scaleArray, dst, scal, nelem, queue);
}
}
