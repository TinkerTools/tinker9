#pragma once
#include "mdprec.h"
#include "time_scale.h"
#include "tool/rc_man.h"


namespace tinker {
void mdrest(int istep);
void mdrest_acc(int istep);


void md_data(rc_op op);


//====================================================================//


void propagate(int nsteps, time_prec dt_ps);


void integrate_data(rc_op);
}


namespace tinker {
void velocity_verlet(int istep, time_prec dt_ps);


// Langevin Piston barostat
// Feller et al. 1995, https://doi.org/10.1063/1.470648
void lpiston_npt(int istep, time_prec dt_ps);
// old xyz
extern pos_prec* leapfrog_x;
extern pos_prec* leapfrog_y;
extern pos_prec* leapfrog_z;
// halftime velocity
extern vel_prec* leapfrog_vx;
extern vel_prec* leapfrog_vy;
extern vel_prec* leapfrog_vz;
// old halftime velocity
extern vel_prec* leapfrog_vxold;
extern vel_prec* leapfrog_vyold;
extern vel_prec* leapfrog_vzold;
extern double hdot_lp;     // box length (h) velocity
extern double hmass_lp;    // h mass
extern double pnhv_lp;     // thermostat velocity
extern double pnhv_pre_lp; // old thermostat velocity
extern double pnhm_lp;     // thermostat mass
extern double pnhf_lp;     // thermostat force
extern double pnh_lp;      // thermostat
void langevin_piston(time_prec dt, virial_prec press);


extern grad_prec *gx1, *gy1, *gz1;
extern grad_prec *gx2, *gy2, *gz2;
constexpr unsigned RESPA_FAST = 1; // 2**0, fast group shall be 0.
constexpr unsigned RESPA_SLOW = 2; // 2**1, slow group shall be 1.
const TimeScaleConfig& respa_tsconfig();
void respa_fast_slow(int istep, time_prec dt_ps);
}
