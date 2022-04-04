#pragma once
#include "ff/energy.h"

namespace tinker {
// Langevin Piston barostat (Leap Frog)
void lf_lpiston_npt(int istep, time_prec dt_ps);
void lf_langevin_piston(time_prec dt, virial_prec press);
extern energy_prec eksum_old; // Kinetic energy at n-1/2.
extern energy_prec eksum_mid; // Kinetic energy at n+1/2.
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

void kinetic_leapfrog(T_prec& temp);

void shake2(time_prec dt, const vel_prec* vxold, const vel_prec* vyold, const vel_prec* vzold,
   const vel_prec* vxnew, const vel_prec* vynew, const vel_prec* vznew, const pos_prec* xold,
   const pos_prec* yold, const pos_prec* zold);
void shake2_acc(time_prec dt, const vel_prec* vxold, const vel_prec* vyold, const vel_prec* vzold,
   const vel_prec* vxnew, const vel_prec* vynew, const vel_prec* vznew, const pos_prec* xold,
   const pos_prec* yold, const pos_prec* zold);

//====================================================================//

void swap_velocity(vel_prec* vxnew, vel_prec* vynew, vel_prec* vznew, vel_prec* vxold,
   vel_prec* vyold, vel_prec* vzold);
void swap_velocity_acc(vel_prec* vxnew, vel_prec* vynew, vel_prec* vznew, vel_prec* vxold,
   vel_prec* vyold, vel_prec* vzold);

void propagate_pos_lp(time_prec dt, pos_prec* x_lp, pos_prec* y_lp, pos_prec* z_lp,
   const vel_prec* vx_lp, const vel_prec* vy, const vel_prec* vz, const pos_prec* xold_lp,
   const pos_prec* yold_lp, const pos_prec* zold_lp, double scale);
void propagate_pos_lp_acc(time_prec dt, pos_prec* x_lp, pos_prec* y_lp, pos_prec* z_lp,
   const vel_prec* vx_lp, const vel_prec* vy_lp, const vel_prec* vz_lp, const pos_prec* xold_lp,
   const pos_prec* yold_lp, const pos_prec* zold_lp, double scale);

void propagate_pos_lp2(time_prec dt, const pos_prec* x_lp, const pos_prec* y_lp,
   const pos_prec* z_lp, pos_prec* xold_lp, pos_prec* yold_lp, pos_prec* zold_lp, double scale);
void propagate_pos_lp2_acc(time_prec dt, const pos_prec* x_lp, const pos_prec* y_lp,
   const pos_prec* z_lp, pos_prec* xold_lp, pos_prec* yold_lp, pos_prec* zold_lp, double scale);

void propagate_pos_lf(time_prec dt, pos_prec* qx, pos_prec* qy, pos_prec* qz, const pos_prec* qxold,
   const pos_prec* qyold, const pos_prec* qzold, const vel_prec* vlx, const vel_prec* vly,
   const vel_prec* vlz);
void propagate_pos_lf_acc(time_prec dt, pos_prec* qx, pos_prec* qy, pos_prec* qz,
   const pos_prec* qxold, const pos_prec* qyold, const pos_prec* qzold, const vel_prec* vlx,
   const vel_prec* vly, const vel_prec* vlz);

void propagate_velocity_lp(vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp,
   const vel_prec* vxnew_lp, const vel_prec* vynew_lp, const vel_prec* vznew_lp,
   const vel_prec* vxold_lp, const vel_prec* vyold_lp, const vel_prec* vzold_lp, const double scale,
   energy_prec& eksum_new, energy_prec& eksum_old);
void propagate_velocity_lp_acc(vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp,
   const vel_prec* vxnew_lp, const vel_prec* vynew_lp, const vel_prec* vznew_lp,
   const vel_prec* vxold_lp, const vel_prec* vyold_lp, const vel_prec* vzold_lp, const double scale,
   energy_prec& eksum_new, energy_prec& eksum_old);

void propagate_velocity_lp2(time_prec dt, vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp,
   const pos_prec* x_lp, const pos_prec* y_lp, const pos_prec* z_lp, const pos_prec* xold_lp,
   const pos_prec* yold_lp, const pos_prec* zold_lp);
void propagate_velocity_lp2_acc(time_prec dt, vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp,
   const pos_prec* x_lp, const pos_prec* y_lp, const pos_prec* z_lp, const pos_prec* xold_lp,
   const pos_prec* yold_lp, const pos_prec* zold_lp);

void propagate_velocity_lp3(vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp,
   const vel_prec* vxnew_lp, const vel_prec* vynew_lp, const vel_prec* vznew_lp,
   const vel_prec* vxold_lp, const vel_prec* vyold_lp, const vel_prec* vzold_lp,
   energy_prec& eksum_new);
void propagate_velocity_lp3_acc(vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp,
   const vel_prec* vxnew_lp, const vel_prec* vynew_lp, const vel_prec* vznew_lp,
   const vel_prec* vxold_lp, const vel_prec* vyold_lp, const vel_prec* vzold_lp,
   energy_prec& eksum_new);
}
