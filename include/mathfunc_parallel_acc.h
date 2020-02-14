#pragma once
#include "dmflag.h"
#include "macro.h"
#include <cstddef>


TINKER_NAMESPACE_BEGIN
namespace platform {
namespace acc {
template <class T>
T reduce_sum(const T* gpu_a, size_t nelem, DMFlag flag);


template <class HT, size_t HN, class DPTR>
void reduce_sum2(HT (&h_ans)[HN], DPTR v, size_t nelem, DMFlag flag);


template <class T>
T reduce_logic_or(const T* a, size_t nelem, DMFlag flag);


template <class T>
T dotprod(const T* a, const T* b, size_t nelem, DMFlag flag);


template <class T>
void dotprod(T* ans, const T* a, const T* b, int nelem, DMFlag flag);


template <class T>
void scale_array(T* dst, T scal, size_t nelem, DMFlag flag);
}
}
TINKER_NAMESPACE_END
