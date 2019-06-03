#ifndef TINKER_GPU_DECL_POLGRP_H_
#define TINKER_GPU_DECL_POLGRP_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
struct polargroup_st {
  static const int maxp11 = 128;
  static const int maxp12 = 128; // 120 -> 128
  static const int maxp13 = 128; // 120 -> 128
  static const int maxp14 = 128; // 120 -> 128

  int *np11, *np12, *np13, *np14;
  int (*ip11)[maxp11];
  int (*ip12)[maxp12];
  int (*ip13)[maxp13];
  int (*ip14)[maxp14];
};

extern polargroup_st polargroup_obj_;
extern polargroup_st* polargroup;

void polargroup_data(int op);
}
TINKER_NAMESPACE_END

#endif