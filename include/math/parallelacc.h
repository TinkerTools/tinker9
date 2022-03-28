#pragma once
#include <cstddef>

namespace tinker {
template <class T>
T reduceSum_acc(const T* gpu_a, size_t nelem, int queue);

template <class HT, size_t HN, class DPTR>
void reduceSum2_acc(HT (&h_ans)[HN], DPTR v, size_t nelem, int queue);

template <class T>
void reduceSumOnDevice_acc(T*, const T*, size_t, int);

template <class HT, size_t HN, class DPTR>
void reduceSum2OnDevice_acc(HT (&)[HN], DPTR, size_t, int);

template <class T>
T dotProd_acc(const T* a, const T* b, size_t nelem, int queue);

template <class T>
void dotProd_acc(T* ans, const T* a, const T* b, size_t nelem, int queue);

template <class T>
void scaleArray_acc(T* dst, T scal, size_t nelem, int queue);
}
