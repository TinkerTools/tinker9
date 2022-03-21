#pragma once
#include "macro.h"

namespace tinker {
/// \ingroup acc_syntax
struct alignas(8) int2
{
   int x, y;
};

/// \ingroup acc_syntax
struct int3
{
   int x, y, z;
};

/// \ingroup acc_syntax
struct alignas(16) int4
{
   int x, y, z, w;
};

/// \ingroup acc_syntax
#pragma acc routine seq
inline int2 make_int2(int x, int y)
{
   int2 f{x, y};
   return f;
}

/// \ingroup acc_syntax
#pragma acc routine seq
inline int3 make_int3(int x, int y, int z)
{
   int3 f{x, y, z};
   return f;
}

/// \ingroup acc_syntax
#pragma acc routine seq
inline int4 make_int4(int x, int y, int z, int w)
{
   int4 f{x, y, z, w};
   return f;
}

/// \ingroup acc_syntax
struct alignas(8) float2
{
   float x, y;
};

/// \ingroup acc_syntax
struct float3
{
   float x, y, z;
};

/// \ingroup acc_syntax
struct alignas(16) float4
{
   float x, y, z, w;
};

/// \ingroup acc_syntax
#pragma acc routine seq
inline float2 make_float2(float x, float y)
{
   float2 f{x, y};
   return f;
}

/// \ingroup acc_syntax
#pragma acc routine seq
inline float3 make_float3(float x, float y, float z)
{
   float3 f{x, y, z};
   return f;
}

/// \ingroup acc_syntax
#pragma acc routine seq
inline float4 make_float4(float x, float y, float z, float w)
{
   float4 f{x, y, z, w};
   return f;
}

/// \ingroup acc_syntax
struct alignas(16) double2
{
   double x, y;
};

/// \ingroup acc_syntax
struct double3
{
   double x, y, z;
};

/// \ingroup acc_syntax
struct alignas(16) double4
{
   double x, y, z, w;
};

/// \ingroup acc_syntax
#pragma acc routine seq
inline double2 make_double2(double x, double y)
{
   double2 f{x, y};
   return f;
}

/// \ingroup acc_syntax
#pragma acc routine seq
inline double3 make_double3(double x, double y, double z)
{
   double3 f{x, y, z};
   return f;
}

/// \ingroup acc_syntax
#pragma acc routine seq
inline double4 make_double4(double x, double y, double z, double w)
{
   double4 f{x, y, z, w};
   return f;
}
}
