#include "md.h"

TINKER_NAMESPACE_BEGIN
size_t estimate_ngangs(int nelem, size_t elem_bytes, size_t max_MB) {
  size_t max_bytes = max_MB * 1024 * 1024;
  size_t denom = nelem * elem_bytes;
  size_t ans = max_bytes / denom;
  ans = nelem ? ans : nelem;

  if (ans >= 0x800)
    return 0x800; // 2048
  else if (ans >= 0x400)
    return 0x400; // 1024
  else if (ans >= 0x200)
    return 0x200; // 512
  else if (ans >= 0x100)
    return 0x100; // 256
  else if (ans >= 0x080)
    return 0x080; // 128
  else
    return 0x040; // 64
}
TINKER_NAMESPACE_END

#if defined(TINKER_DOUBLE_PRECISION)
#  include "md_egv_d.hh"
#elif defined(TINKER_SINGLE_PRECISION)
#  include "md_egv_s.hh"
#else
static_assert(false, "");
#endif
