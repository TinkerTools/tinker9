#pragma once
#include <cstddef>

namespace tinker {
inline void host_zero(int& v)
{
   v = 0;
}

inline void host_zero(float& v)
{
   v = 0;
}

inline void host_zero(double& v)
{
   v = 0;
}

inline void host_zero(unsigned long long& v)
{
   v = 0;
}

template <class T>
void host_zero(T*& ptr)
{
   ptr = nullptr;
}

template <class T, size_t N>
void host_zero(T (&v)[N])
{
   for (size_t i = 0; i < N; ++i)
      host_zero(v[i]);
}

template <class T, class... Ts>
void host_zero(T& v, Ts&... vs)
{
   host_zero(v);
   host_zero(vs...);
}
}
