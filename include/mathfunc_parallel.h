#pragma once

#include "deduce_ptr.h"

TINKER_NAMESPACE_BEGIN
/// Functions running in parallel.
namespace parallel {
/**
 * \ingroup math
 * \brief Sum all the elements of an array.
 * \return
 * The sum.
 */
template <class T>
T reduce_sum(const T* gpu_a, size_t nelem);

template <class HT, size_t HN, class DPTR>
void reduce_sum2(HT (&h_ans)[HN], DPTR v, size_t nelem);

/**
 * \ingroup math
 * \brief N-dimensional dot product.
 * \f[ DotProduct = \sum_i^n a_i \cdot b_i \f]
 *
 * \param[in] a
 * Device pointer to array a.
 * \param[in] b
 * Device pointer to array b.
 * \param[in] nelem
 * Number of elements in each array.
 *
 * @return
 * The dot product to the host thread.
 */
float dotprod(const float* a, const float* b, size_t nelem);
/// \ingroup math
double dotprod(const double* a, const double* b, size_t nelem);

/**
 * \ingroup math
 * \brief Multiply all of the elements in an array by a scalar.
 * \f[ a_i = c \cdot a_i, 1 \leq i \leq n \f]
 *
 * \param[in,out] dst
 * Device pointer to the array.
 * \param[in] scal
 * The multiplier, a scalar.
 * \param[in] nelem
 * Number of elements in the array.
 */
void scale_array(float* dst, float scal, size_t nelem);
/// \ingroup math
void scale_array(double* dst, double scal, size_t nelem);
}
TINKER_NAMESPACE_END
