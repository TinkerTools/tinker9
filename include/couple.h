#pragma once
#include "tool/rcman.h"

namespace tinker {
const int couple_maxn12 = 8;
extern int (*couple_i12)[couple_maxn12];
extern int* couple_n12;

void couple_data(RcOp);
}
