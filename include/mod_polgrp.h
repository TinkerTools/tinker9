#ifndef TINKER_MOD_POLGRP_H_
#define TINKER_MOD_POLGRP_H_

#include "util_cxx.h"

TINKER_NAMESPACE_BEGIN
struct PolarizationGroupConnectivityLists {
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
typedef PolarizationGroupConnectivityLists polargroup_t;

TINKER_EXTERN polargroup_t polargroup_obj_;
TINKER_EXTERN polargroup_t* polargroup;

TINKER_NAMESPACE_END

#endif
