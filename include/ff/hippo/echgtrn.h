#pragma once
#include "ff/elec.h"
#include "glob/mplpot.h"
#include "glob/mpole.h"
#include "mod/elechippo.h"
#include "tool/rcman.h"

namespace tinker {
void echgtrn_data(RcOp);
void echgtrn(int vers);
void echgtrn_cu(int);
void echgtrn_acc(int);
}
