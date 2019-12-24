#pragma once
#include "realn_def.h"
#include "seq_def.h"


TINKER_NAMESPACE_BEGIN
#if TINKER_DOUBLE_PRECISION
using real2 = double2;
using real3 = double3;
using real4 = double4;
#   define make_real2(x, y) make_double2((x), (y))
#   define make_real3(x, y, z) make_double3((x), (y), (z))
#   define make_real4(x, y, z, w) make_double4((x), (y), (z), (w))
#endif


#if TINKER_SINGLE_PRECISION
using real2 = float2;
using real3 = float3;
using real4 = float4;
#   define make_real2(x, y) make_float2((x), (y))
#   define make_real3(x, y, z) make_float3((x), (y), (z))
#   define make_real4(x, y, z, w) make_float4((x), (y), (z), (w))

#endif


// + - *
ROUTINE_SEQ
inline real3 operator+(real3 a, real3 b)
{
   return make_real3(a.x + b.x, a.y + b.y, a.z + b.z);
}


ROUTINE_SEQ
inline real3 operator-(real3 a, real3 b)
{
   return make_real3(a.x - b.x, a.y - b.y, a.z - b.z);
}


ROUTINE_SEQ
inline real3 operator*(real a, real3 b)
{
   return make_real3(a * b.x, a * b.y, a * b.z);
}


// += -=
ROUTINE_SEQ
inline void operator+=(real3& a, real3 b)
{
   a.x += b.x;
   a.y += b.y;
   a.z += b.z;
}


ROUTINE_SEQ
inline void operator-=(real3& a, real3 b)
{
   a.x -= b.x;
   a.y -= b.y;
   a.z -= b.z;
}


// dot
ROUTINE_SEQ
inline real dot3(real3 a, real3 b)
{
   return a.x * b.x + a.y * b.y + a.z * b.z;
}


ROUTINE_SEQ
inline real dot3(real ax, real ay, real az, real bx, real by, real bz)
{
   return ax * bx + ay * by + az * bz;
}
TINKER_NAMESPACE_END
