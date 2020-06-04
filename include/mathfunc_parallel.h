#pragma once
#include "mathfunc_parallel_acc.h"
#include "mathfunc_parallel_cu.h"
#include "platform.h"


namespace tinker {
namespace parallel {
/**
 * \ingroup math
 * \brief Sum over all of the elements of an 1D array.
 *
 * \f[ Sum = \sum_i^n a_i \f]
 * \return The sum.
 */
template <class T>
T reduce_sum(const T* gpu_a, size_t nelem, LPFlag flag)
{
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      return reduce_sum_cu(gpu_a, nelem, flag);
   else
#endif
      return reduce_sum_acc(gpu_a, nelem, flag);
}


/**
 * \ingroup math
 * \brief Sum over all of the elements of a 2D array.
 *
 * Fortran syntax:
 * \f[ Ans(k) = \sum_i^n v(k,i), 1 \le k \le HN \f]
 * C++ syntax:
 * \f[ Ans[k] = \sum_i^n v[i][k], 0 \le k < HN \f]
 */
template <class HT, size_t HN, class DPTR>
void reduce_sum2(HT (&h_ans)[HN], DPTR v, size_t nelem, LPFlag flag)
{
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      reduce_sum2_cu(h_ans, v, nelem, flag);
   else
#endif
      reduce_sum2_acc(h_ans, v, nelem, flag);
}


template <class T>
T reduce_logic_or(const T* a, size_t nelem, LPFlag flag)
{
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      return reduce_logic_or_cu(a, nelem, flag);
   else
#endif
      return reduce_logic_or_acc(a, nelem, flag);
}


/**
 * \ingroup math
 * \brief Dot product of two linear arrays.
 *
 * \f[ DotProduct = \sum_i^n a_i \cdot b_i \f]
 * \return The dot product to the host thread.
 */
template <class T>
T dotprod(const T* a, const T* b, size_t nelem, LPFlag flag)
{
   return dotprod_acc(a, b, nelem, flag);
}


/**
 * \ingroup math
 * \brief Dot product of two linear arrays.
 */
template <class T>
void dotprod(T* ans, const T* a, const T* b, int nelem, LPFlag flag)
{
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      dotprod_cu(ans, a, b, nelem, flag);
   else
#endif
      dotprod_acc(ans, a, b, nelem, flag);
}


/**
 * \ingroup math
 * \brief Multiply all of the elements in an 1D array by a scalar.
 *
 * \f[ a_i = c \cdot a_i \f]
 */
template <class T>
void scale_array(T* dst, T scal, size_t nelem, LPFlag flag)
{
   return scale_array_acc(dst, scal, nelem, flag);
}
}
}
