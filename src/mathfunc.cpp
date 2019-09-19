#include "mathfunc.h"
#include <cassert>

TINKER_NAMESPACE_BEGIN
bool is_pow2(size_t x) { return (x != 0) && ((x & (x - 1)) == 0); }

size_t pow2_le(size_t x) {
  assert(val > 0);
  return 1ull << (sizeof(size_t) * 8 - 1 - __builtin_clzll(x));
}

size_t pow2_ge(size_t x) {
  if (x <= 1)
    return 1;
  return 1ull << (sizeof(size_t) * 8 - __builtin_clzll(x - 1));
}
TINKER_NAMESPACE_END
