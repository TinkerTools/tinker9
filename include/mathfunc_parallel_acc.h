#pragma once
#include "macro.h"
#include <cstddef>


namespace tinker {
template <class T>
T reduce_sum_acc(const T* gpu_a, size_t nelem, int queue);


template <class HT, size_t HN, class DPTR>
void reduce_sum2_acc(HT (&h_ans)[HN], DPTR v, size_t nelem, int queue);


/**
 * \ingroup math
 * \see reduce_sum_on_device
 */
template <class T>
void reduce_sum_on_device_acc(T* dp_ans, T& ans, const T* a, size_t nelem,
                              int queue);


template <class T>
T dotprod_acc(const T* a, const T* b, size_t nelem, int queue);


template <class T>
void dotprod_acc(T* ans, const T* a, const T* b, size_t nelem, int queue);


template <class T>
void scale_array_acc(T* dst, T scal, size_t nelem, int queue);
}
