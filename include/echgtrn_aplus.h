#pragma once
#include "elec.h"
#include "mod.chgtrn.h"
#include "mod.ctrpot.h"
#include "mod.mplpot.h"
#include "mod.mpole.h"
#include "tool/rc_man.h"


namespace tinker {
void echgtrn_aplus_data(rc_op op);
void echgtrn_aplus(int vers);
void echgtrn_aplus_cu(int);
void echgtrn_aplus_acc(int);
}
