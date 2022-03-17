#pragma once
#include <cstddef>

namespace tinker {
//====================================================================//
// host

inline void zeroOnHost(int& v)
{
   v = 0;
}

inline void zeroOnHost(float& v)
{
   v = 0;
}

inline void zeroOnHost(double& v)
{
   v = 0;
}

inline void zeroOnHost(unsigned long long& v)
{
   v = 0;
}

template <class T>
void zeroOnHost(T*& ptr)
{
   ptr = nullptr;
}

template <class T, size_t N>
void zeroOnHost(T (&v)[N])
{
   for (size_t i = 0; i < N; ++i)
      zeroOnHost(v[i]);
}

template <class T, class... Ts>
void zeroOnHost(T& v, Ts&... vs)
{
   zeroOnHost(v);
   zeroOnHost(vs...);
}

//====================================================================//
// device

template <class T>
void zeroOnDevice3Async_acc(int nelem, T* a1, T* a2, T* a3);

template <class T>
void zeroOnDevice3Async(int nelem, T* a1, T* a2, T* a3)
{
   zeroOnDevice3Async_acc(nelem, a1, a2, a3);
}

template <class T, int N>
void zeroOnDevice3Async_acc(int nelem, T (*a1)[N], T (*a2)[N], T (*a3)[N]);

template <class T, int N>
void zeroOnDevice3Async(int nelem, T (*a1)[N], T (*a2)[N], T (*a3)[N])
{
   zeroOnDevice3Async_acc<T, N>(nelem, a1, a2, a3);
}

template <class T>
void zeroOnDevice9Async_acc(
   int nelem, T* a1, T* a2, T* a3, T* a4, T* a5, T* a6, T* a7, T* a8, T* a9);

template <class T>
void zeroOnDevice9Async(int nelem, T* a1, T* a2, T* a3, T* a4, T* a5, T* a6, T* a7, T* a8, T* a9)
{
   zeroOnDevice9Async_acc(nelem, a1, a2, a3, a4, a5, a6, a7, a8, a9);
}
}
