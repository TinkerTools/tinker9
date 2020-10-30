#pragma once
#include "macro.h"
#include "mdprec.h"
#include "tool/rc_man.h"


namespace tinker {
void accmanaged_data(rc_op);


namespace detail {
extern energy_prec exx, eyy, ezz, exy, eyz, ezx; // kinetic
extern vel_prec vtot1, vtot2, vtot3;             // mdrest
}
}
