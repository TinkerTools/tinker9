#ifndef TINKER_COUPLE_H_
#define TINKER_COUPLE_H_

#include "gen_unit.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
/// @brief
/// atom neighbor connectivity lists
struct Couple {
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

typedef GenericUnit<Couple, 1> CoupleUnit;
TINKER_EXTERN CoupleUnit couple_unit;

void couple_data(rc_op);
TINKER_NAMESPACE_END

#endif
