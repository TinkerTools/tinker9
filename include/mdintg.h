#pragma once
#include "mdprec.h"
#include "time_scale.h"
#include "tool/rc_man.h"


namespace tinker {
void mdrest(int istep);
void mdrest_acc(int istep);
void mdrest_remove_pbc_momentum_cu(bool copyout, vel_prec& vtot1,
                                   vel_prec& vtot2, vel_prec& vtot3);


void md_data(rc_op op);


//====================================================================//


void propagate(int nsteps, time_prec dt_ps);


void integrate_data(rc_op);
}


namespace tinker {
void velocity_verlet(int istep, time_prec dt_ps);


extern grad_prec *gx1, *gy1, *gz1;
extern grad_prec *gx2, *gy2, *gz2;
constexpr unsigned RESPA_FAST = 1; // 2**0, fast group shall be 0.
constexpr unsigned RESPA_SLOW = 2; // 2**1, slow group shall be 1.
const TimeScaleConfig& respa_tsconfig();
void respa_fast_slow(int istep, time_prec dt_ps);


void vv_lpiston_npt(int istep, time_prec dt_ps);
}
