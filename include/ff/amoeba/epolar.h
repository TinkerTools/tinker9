#pragma once
#include "ff/elec.h"
#include "tool/rcman.h"

namespace tinker {
void epolarData(RcOp);
void epolar(int vers);
// see also subroutine epolar0e in epolar.f
void epolar0DotProd(const real (*uind)[3], const real (*udirp)[3]);
}
