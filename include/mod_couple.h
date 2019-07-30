#ifndef TINKER_MOD_COUPLE_H_
#define TINKER_MOD_COUPLE_H_

#include "util_cxx.h"

TINKER_NAMESPACE_BEGIN
struct couple_t {
  static const int maxn12 = 8;
  static const int maxn13 = 32;  // 24 -> 32
  static const int maxn14 = 96;  // 72 -> 96
  static const int maxn15 = 224; // 216 -> 224

  int *n12, *n13, *n14, *n15;
  int (*i12)[maxn12];
  int (*i13)[maxn13];
  int (*i14)[maxn14];
  int (*i15)[maxn15];
};

TINKER_EXTERN couple_t coupl_obj_;
TINKER_EXTERN couple_t* coupl;
TINKER_NAMESPACE_END

#endif
