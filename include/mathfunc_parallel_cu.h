#pragma once
#include "macro.h"
#include <cstddef>


namespace tinker {
template <class T>
T reduce_sum_cu(const T* a, size_t nelem, int queue);


template <class HT, size_t HN, class DPTR>
void reduce_sum2_cu(HT (&h_ans)[HN], DPTR v, size_t nelem, int queue);


/**
 * \ingroup math
 * \see parallel::reduce_sum_on_device
 */
template <class T>
void reduce_sum_on_device_cu(T* dp_ans, T& ans, const T* a, size_t nelem,
                             int queue);


template <class T>
void dotprod_cu(T* ans, const T* a, const T* b, size_t nelem, int queue);
}
