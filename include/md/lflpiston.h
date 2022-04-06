#pragma once
#include "ff/precision.h"

namespace tinker {
void lf_lpiston_npt(int istep, time_prec dt_ps); // Langevin Piston barostat (Leap Frog)
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
TINKER_EXTERN energy_prec eksum_old; // Kinetic energy at n-1/2.
TINKER_EXTERN energy_prec eksum_mid; // Kinetic energy at n+1/2.
// old xyz
TINKER_EXTERN pos_prec* leapfrog_x;
TINKER_EXTERN pos_prec* leapfrog_y;
TINKER_EXTERN pos_prec* leapfrog_z;
// halftime velocity
TINKER_EXTERN vel_prec* leapfrog_vx;
TINKER_EXTERN vel_prec* leapfrog_vy;
TINKER_EXTERN vel_prec* leapfrog_vz;
// old halftime velocity
TINKER_EXTERN vel_prec* leapfrog_vxold;
TINKER_EXTERN vel_prec* leapfrog_vyold;
TINKER_EXTERN vel_prec* leapfrog_vzold;
TINKER_EXTERN double hdot_lp;     // box length (h) velocity
TINKER_EXTERN double hmass_lp;    // h mass
TINKER_EXTERN double pnhv_lp;     // thermostat velocity
TINKER_EXTERN double pnhv_pre_lp; // old thermostat velocity
TINKER_EXTERN double pnhm_lp;     // thermostat mass
TINKER_EXTERN double pnhf_lp;     // thermostat force
TINKER_EXTERN double pnh_lp;      // thermostat
}
