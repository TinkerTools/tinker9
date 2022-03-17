#pragma once
#include "elec.h"
#include "mod.chgtrn.h"
#include "mod.ctrpot.h"
#include "mod.mplpot.h"
#include "mod.mpole.h"
#include "tool/rcman.h"

namespace tinker {
void echgtrn_data(RcOp);
void echgtrn(int vers);
void echgtrn_cu(int);
void echgtrn_acc(int);
}
