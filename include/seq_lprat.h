#pragma once
#include "add.h"
#include "seq_def.h"


namespace tinker {
#pragma acc routine seq
SEQ_CUDA
double lprat_e1(double dt2, double mfrac, double vbar)
{
   return exp(-dt2 * mfrac * vbar);
}


#pragma acc routine seq
SEQ_CUDA
double lprat_e2(double dt2, double mfrac, double vbar, double al, double vnh)
{
   double t1 = 1.0 / (1.0 + mfrac * al * vbar * dt2);
   double t2 = dt2 * (vnh + al * vbar * mfrac);
   return t1 * exp(t2);
}
}
