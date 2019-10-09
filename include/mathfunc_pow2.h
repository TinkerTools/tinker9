#pragma once
#include "macro.h"
#include <cstring>


TINKER_NAMESPACE_BEGIN
/// \ingroup math
/// \return
/// True if and only if \c val is a power of 2.
constexpr bool is_pow2(size_t val)
{
   return (val != 0) && ((val & (val - 1)) == 0);
}


/// \ingroup math
/// \return
/// The base-2 logarithm of positive integer \c n.
/// \param val
/// Must be greater than 0.
constexpr size_t int_log2(size_t val)
{
   return val < 2 ? 0 : 1 + int_log2(val / 2);
}


/// \ingroup math
/// \return
/// A power of 2 that is less than or equal to \c val.
/// \param val
/// Must be greater than 0.
constexpr size_t pow2_le(size_t val)
{
   return 1ull << int_log2(val);
}


/// \ingroup math
/// \return
/// A power of 2 that is greater than or equal to \c val.
constexpr size_t pow2_ge(size_t val)
{
   return val <= 1 ? 1 : (1ull << (1 + int_log2(val - 1)));
}
TINKER_NAMESPACE_END
