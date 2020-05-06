#pragma once
#include "dmflag.h"
#include <cstddef>


namespace tinker {
template <class T>
T reduce_sum_cu(const T* a, size_t nelem, DMFlag flag);


template <class HT, size_t HN, class DPTR>
void reduce_sum2_cu(HT (&h_ans)[HN], DPTR v, size_t nelem, DMFlag flag);


template <class T>
T reduce_logic_or_cu(const T* a, size_t nelem, DMFlag flag);


template <class T>
void dotprod_cu(T* ans, const T* a, const T* b, int nelem, DMFlag flag);
}
