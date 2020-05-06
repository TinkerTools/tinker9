#pragma once
#include "dmflag.h"
#include "macro.h"
#include <cstddef>


namespace tinker {
template <class T>
T reduce_sum_acc(const T* gpu_a, size_t nelem, DMFlag flag);


template <class HT, size_t HN, class DPTR>
void reduce_sum2_acc(HT (&h_ans)[HN], DPTR v, size_t nelem, DMFlag flag);


template <class T>
T reduce_logic_or_acc(const T* a, size_t nelem, DMFlag flag);


template <class T>
T dotprod_acc(const T* a, const T* b, size_t nelem, DMFlag flag);


template <class T>
void dotprod_acc(T* ans, const T* a, const T* b, int nelem, DMFlag flag);


template <class T>
void scale_array_acc(T* dst, T scal, size_t nelem, DMFlag flag);
}
