#ifndef TINKER_GPU_DECL_COUPLE_H_
#define TINKER_GPU_DECL_COUPLE_H_

#include "decl_basic.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
struct couple_st {
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

extern couple_st couple_obj_;
extern couple_st* couple;

void couple_data(rc_t rc);
}
TINKER_NAMESPACE_END

#endif
