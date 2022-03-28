#pragma once
#include "macro.h"
#include "platform.h"

#include "math/parallelacc.h"
#include "math/parallelcu.h"

namespace tinker {
/// \ingroup math_parallel
/// \brief Sum over all of the elements of an 1D array.
///
/// \f[ Sum = \sum_i^n a_i \f]
/// \return The sum.
template <class T>
T reduceSum(const T* gpu_a, size_t nelem, int queue)
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      return reduceSum_cu(gpu_a, nelem, queue);
   else
#endif
      return reduceSum_acc(gpu_a, nelem, queue);
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
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      reduceSum2_cu(h_ans, v, nelem, queue);
   else
#endif
      reduceSum2_acc(h_ans, v, nelem, queue);
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
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      reduceSumOnDevice_cu(dp_ans, a, nelem, queue);
   else
#endif
      reduceSumOnDevice_acc(dp_ans, a, nelem, queue);
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
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      reduceSum2OnDevice_cu(dref, v, nelem, queue);
   else
#endif
      reduceSum2OnDevice_acc(dref, v, nelem, queue);
}

/// \ingroup math_parallel
/// \brief Dot product of two linear arrays.
///
/// \f[ DotProduct = \sum_i^n a_i \cdot b_i \f]
/// \return The dot product to the host thread.
template <class T>
T dotProd(const T* a, const T* b, size_t nelem, int queue)
{
   return dotProd_acc(a, b, nelem, queue);
}

/// \ingroup math_parallel
/// \brief Dot product of two linear arrays.
template <class T>
void dotProd(T* ans, const T* a, const T* b, size_t nelem, int queue)
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      dotProd_cu(ans, a, b, nelem, queue);
   else
#endif
      dotProd_acc(ans, a, b, nelem, queue);
}

/// \ingroup math_parallel
/// \brief Multiply all of the elements in an 1D array by a scalar.
///
/// \f[ a_i = c \cdot a_i \f]
template <class T>
void scaleArray(T* dst, T scal, size_t nelem, int queue)
{
   return scaleArray_acc(dst, scal, nelem, queue);
}
}
