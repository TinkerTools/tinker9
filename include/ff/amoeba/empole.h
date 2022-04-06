#pragma once
#include "ff/precision.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup mpole
/// \{
void empoleData(RcOp);
void empole(int vers);
void empoleEwaldRecip(int vers);
void torque(int vers, grad_prec* dx, grad_prec* dy, grad_prec* dz);
void mpoleInit(int vers);
/// \}
}
