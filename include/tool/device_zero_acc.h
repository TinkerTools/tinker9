#pragma once
#include "macro.h"


namespace tinker {
template <class T>
void zero3_async_acc(int n, T* a1, T* a2, T* a3);
template <class T, int N>
void zero3_async_acc(int nelem, T (*a1)[N], T (*a2)[N], T (*a3)[N]);


template <class T>
void zero9_async_acc(int nelem, T* a1, T* a2, T* a3, T* a4, T* a5, T* a6, T* a7,
                     T* a8, T* a9);
}
