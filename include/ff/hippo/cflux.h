#pragma once
#include "ff/energybuffer.h"
#include "tool/rcman.h"

namespace tinker {
void cfluxData(RcOp);

void alterchg();
void cfluxZeroPot();
void dcflux(int vers, grad_prec* gx, grad_prec* gy, grad_prec* gz, VirialBuffer vir);
}
