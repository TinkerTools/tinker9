#pragma once
#include "ff/precision.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup hippopolar
void epolarChgpenData(RcOp);
/// \ingroup hippopolar
void epolarChgpen(int vers);
/// \ingroup hippopolar
void epolarChgpenEwaldRecipSelf(int vers, int use_cf);
}
