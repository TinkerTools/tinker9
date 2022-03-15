#pragma once
#include "nhc.h"

namespace tinker {
// Martyna et al. 1996, https://doi.org/10.1080/00268979600100761
void nhc_npt(int istep, time_prec dt_ps);
// section 4.4: iLnhc + iLp
void hoover(time_prec dt, virial_prec press);

extern double vbar;         // epsilon velocity
extern double qbar;         // epsilon mass
extern double gbar;         // epsilon force
extern double vnh[maxnose]; // ksi velocity
extern double qnh[maxnose]; // ksi mass
extern double gnh[maxnose]; // ksi force
}
