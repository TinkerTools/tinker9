#pragma once
#include "mdprec.h"

namespace tinker {
void __debug_norm_propagate_pos_acc(
   pos_prec poseps, time_prec dt, const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);
}
