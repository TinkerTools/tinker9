#pragma once
#include "ff/time_scale.h"

namespace tinker {
void mdrest(int istep);
void mdData(RcOp);
void mdPropagate(int nsteps, time_prec dt_ps);
void mdIntegrateData(RcOp);

constexpr unsigned RESPA_FAST = 1; // 2**0, fast group shall be 0.
constexpr unsigned RESPA_SLOW = 2; // 2**1, slow group shall be 1.
const TimeScaleConfig& mdRespaTsconfig();

void mdsaveAsync(int istep, time_prec dt);
void mdsaveSynchronize();
void mdsaveData(RcOp);
}
