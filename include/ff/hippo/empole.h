#pragma once
#include "tool/rcman.h"

namespace tinker {
/// \ingroup hippompole
void empoleChgpenData(RcOp);
/// \ingroup hippompole
void empoleChgpen(int vers);
/// \ingroup hippompole
void empoleChgpenEwaldRecip(int vers, int useCF);
}

namespace tinker {
void empoleAplusEwald(int vers, int useCF);
void empoleAplusNonEwald(int vers, int useCF);
}
