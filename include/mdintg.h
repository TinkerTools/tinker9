#pragma once
#include "mdprec.h"
#include "time_scale.h"
#include "tool/rc_man.h"

namespace tinker {
void mdrest(int istep);
void md_data(rc_op op);
void propagate(int nsteps, time_prec dt_ps);
void integrate_data(rc_op);

extern grad_prec *gx1, *gy1, *gz1;
extern grad_prec *gx2, *gy2, *gz2;
constexpr unsigned RESPA_FAST = 1; // 2**0, fast group shall be 0.
constexpr unsigned RESPA_SLOW = 2; // 2**1, slow group shall be 1.
const TimeScaleConfig& respa_tsconfig();
}
