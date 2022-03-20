#pragma once
#include "elec.h"
#include "glob/chgtrn.h"
#include "glob/ctrpot.h"
#include "glob/mplpot.h"
#include "glob/mpole.h"
#include "tool/rcman.h"

namespace tinker {
void echgtrn_data(RcOp);
void echgtrn(int vers);
void echgtrn_cu(int);
void echgtrn_acc(int);
}
