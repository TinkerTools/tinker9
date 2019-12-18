#pragma once
#include "macro.h"
#include "macro_void_cuda_def.h"


TINKER_NAMESPACE_BEGIN
#ifndef __CUDACC__
#pragma acc routine seq
inline float2 make_float2(float x, float y)
{
   float2 f{.x = x, .y = y};
   return f;
}


#pragma acc routine seq
inline float3 make_float3(float x, float y, float z)
{
   float3 f{.x = x, .y = y, .z = z};
   return f;
}


#pragma acc routine seq
inline float4 make_float4(float x, float y, float z, float w)
{
   float4 f{.x = x, .y = y, .z = z, .w = w};
   return f;
}


#pragma acc routine seq
inline double2 make_double2(double x, double y)
{
   double2 f{.x = x, .y = y};
   return f;
}


#pragma acc routine seq
inline double3 make_double3(double x, double y, double z)
{
   double3 f{.x = x, .y = y, .z = z};
   return f;
}


#pragma acc routine seq
inline double4 make_double4(double x, double y, double z, double w)
{
   double4 f{.x = x, .y = y, .z = z, .w = w};
   return f;
}
#endif


#if TINKER_SINGLE_PRECISION
#   define make_real2(x, y) make_float2((x), (y))
#   define make_real3(x, y, z) make_float3((x), (y), (z))
#   define make_real4(x, y, z, w) make_float4((x), (y), (z), (w))
#endif


#if TINKER_DOUBLE_PRECISION
#   define make_real2(x, y) make_double2((x), (y))
#   define make_real3(x, y, z) make_double3((x), (y), (z))
#   define make_real4(x, y, z, w) make_double4((x), (y), (z), (w))
#endif


// + - *
#pragma acc routine seq
__device__
inline real3 operator+(real3 a, real3 b)
{
   return make_real3(a.x + b.x, a.y + b.y, a.z + b.z);
}


#pragma acc routine seq
__device__
inline real3 operator-(real3 a, real3 b)
{
   return make_real3(a.x - b.x, a.y - b.y, a.z - b.z);
}


#pragma acc routine seq
__device__
inline real3 operator*(real a, real3 b)
{
   return make_real3(a * b.x, a * b.y, a * b.z);
}


// += -=
#pragma acc routine seq
__device__
inline void operator+=(real3& a, real3 b)
{
   a.x += b.x;
   a.y += b.y;
   a.z += b.z;
}


#pragma acc routine seq
__device__
inline void operator-=(real3& a, real3 b)
{
   a.x -= b.x;
   a.y -= b.y;
   a.z -= b.z;
}


// dot
#pragma acc routine seq
__device__
inline real dot3(real3 a, real3 b)
{
   return a.x * b.x + a.y * b.y + a.z * b.z;
}


#pragma acc routine seq
__device__
inline real dot3(real ax, real ay, real az, real bx, real by, real bz)
{
   return ax * bx + ay * by + az * bz;
}
TINKER_NAMESPACE_END
