#pragma once
#include "macro.h"
#include <cstring>


TINKER_NAMESPACE_BEGIN
/// \ingroup math
/// \brief Return true if and only if `val` is a power of 2.
constexpr bool is_pow2(size_t val)
{
   return (val != 0) && ((val & (val - 1)) == 0);
}


/// \ingroup math
/// \brief Integer base 2 power.
constexpr int pow2(int val)
{
   return 1 << val;
}


/// \ingroup math
/// \brief Long long integer base 2 power.
constexpr long long pow2ll(int val)
{
   return 1ll << val;
}


/// \ingroup math
/// \brief The base-2 logarithm of `val`.
/// \param val
/// Must be greater than 0.
/// \note
/// Should return 0 if `val` is 0.
constexpr int int_log2(size_t val)
{
   return val < 2 ? 0 : 1 + int_log2(val >> 1);
}


/// \ingroup math
/// \brief Non-constexpr `int_log2(int)`.
inline int builtin_floor_log2(int val)
{
   if (val == 0)
      return 0;
   int c = __builtin_clz(val);
   return 8 * sizeof(int) - 1 - c;
}


/// \ingroup
/// \brief Non-constexpr `int_log2(long long)`.
inline int builtin_floor_log2(long long val)
{
   if (val == 0)
      return 0;
   int c = __builtin_clzll(val);
   return 8 * sizeof(long long) - 1 - c;
}


/// \ingroup math
/// \brief Non-constexpr log base 2 ceiling.
inline int builtin_ceil_log2(int val)
{
   return builtin_floor_log2(val) + (val & (val - 1) ? 1 : 0);
}


/// \ingroup math
/// \brief Non-constexpr log base 2 ceiling.
inline int builtin_ceil_log2(long long val)
{
   return builtin_floor_log2(val) + (val & (val - 1) ? 1 : 0);
}


/// \ingroup math
/// \brief A power of 2 that is less than or equal to `val`.
/// \param val
/// Must be greater than 0.
constexpr size_t pow2_le(size_t val)
{
   return 1ull << int_log2(val);
}


/// \ingroup math
/// \brief A power of 2 that is greater than or equal to `val`.
constexpr size_t pow2_ge(size_t val)
{
   return val <= 1 ? 1 : (1ull << (1 + int_log2(val - 1)));
}
TINKER_NAMESPACE_END
