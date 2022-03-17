#pragma once

namespace tinker {
template <class T>
void zero3_async_acc(int nelem, T* a1, T* a2, T* a3);

template <class T>
void zero3_async(int nelem, T* a1, T* a2, T* a3)
{
   zero3_async_acc(nelem, a1, a2, a3);
}

template <class T, int N>
void zero3_async_acc(int nelem, T (*a1)[N], T (*a2)[N], T (*a3)[N]);

template <class T, int N>
void zero3_async(int nelem, T (*a1)[N], T (*a2)[N], T (*a3)[N])
{
   zero3_async_acc<T, N>(nelem, a1, a2, a3);
}

template <class T>
void zero9_async_acc(int nelem, T* a1, T* a2, T* a3, T* a4, T* a5, T* a6, T* a7, T* a8, T* a9);

template <class T>
void zero9_async(int nelem, T* a1, T* a2, T* a3, T* a4, T* a5, T* a6, T* a7, T* a8, T* a9)
{
   zero9_async_acc(nelem, a1, a2, a3, a4, a5, a6, a7, a8, a9);
}

}
