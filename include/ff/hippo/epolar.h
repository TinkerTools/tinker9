#pragma once
#include "tool/rcman.h"

namespace tinker {
/// \ingroup hippopolar
void epolarChgpenData(RcOp);
/// \ingroup hippopolar
void epolarChgpen(int vers);
/// \ingroup hippopolar
void epolarChgpenEwaldRecipSelf(int vers, int use_cf);
}

namespace tinker {
void epolarAplusEwald(int vers, int useCF);
void epolarAplusNonEwald(int vers, int useCF);
}
