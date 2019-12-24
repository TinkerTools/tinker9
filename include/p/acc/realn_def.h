#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
struct alignas(8) float2
{
   float x, y;
};


struct float3
{
   float x, y, z;
};


struct alignas(16) float4
{
   float x, y, z, w;
};


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


struct alignas(16) double2
{
   double x, y;
};


struct double3
{
   double x, y, z;
};


struct alignas(16) double4
{
   double x, y, z, w;
};


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
TINKER_NAMESPACE_END
