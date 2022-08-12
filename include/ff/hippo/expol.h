#pragma once
#include "ff/energybuffer.h"
#include "tool/rcman.h"

namespace tinker {
void expolData(RcOp);

void alterpol(real (*polscale)[3][3], real (*polinv)[3][3]);
void dexpol(int vers);
}
