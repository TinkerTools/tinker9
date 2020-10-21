#pragma once
#include "elec.h"
#include "mod.chgtrn.h"
#include "mod.ctrpot.h"
#include "mod.mplpot.h"
#include "mod.mpole.h"
#include "tool/rc_man.h"


namespace tinker {
void echgtrn_data(rc_op op);
void echgtrn(int vers);
void echgtrn_cu(int);
}
