#pragma once
#include "tool/energybuffer.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup cflux
/// \{
void cfluxData(RcOp);

void alterchg();
void cfluxZeroPot();
void dcflux(int vers, grad_prec* gx, grad_prec* gy, grad_prec* gz, VirialBuffer vir);
/// \}
}
