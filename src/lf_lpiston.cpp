#include "lf_lpiston.h"
#include "box.h"
#include "energy.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include "mdpt.h"
#include "random.h"
#include "rattle.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/stodyn.hh>
#include <tinker/detail/units.hh>


namespace tinker {
energy_prec eksum_old;
energy_prec eksum_mid;
pos_prec *leapfrog_x, *leapfrog_y, *leapfrog_z;
vel_prec *leapfrog_vx, *leapfrog_vy, *leapfrog_vz;
vel_prec *leapfrog_vxold, *leapfrog_vyold, *leapfrog_vzold;
double hdot_lp;
double hmass_lp;
double pnhv_lp;
double pnhv_pre_lp;
double pnhm_lp;
double pnhf_lp;
double pnh_lp;


namespace {
double press;
}


void lf_langevin_piston(time_prec dt, virial_prec press)
{
   const double vbox = volbox();
   const int df = mdstuf::nfree;
   hmass_lp =
      bath::taupres; // use these as input for the baro/thermostat masses
   pnhm_lp = bath::tautemp;

   double dpress;

   double xtlabc = lvec1.x;
   double delpr;
   double delpx;
   double hdot_old;
   double hdot_pre;
   double hdot_mid;
   double sf, sh;
   const double gamma_piston = stodyn::friction;
   const double palpha =
      (1.0 - gamma_piston * dt / 2.0) / (1.0 + gamma_piston * dt / 2.0);
   const double pbfact = dt * dt / (1.0 + gamma_piston * dt / 2.0);
   const double kbt = units::gasconst * bath::kelvin;
   const double prfwd =
      sqrt(2.0 * gamma_piston * dt * kbt / (3.0 * hmass_lp)) / dt;
   double h_random_force;
   double pwinv;
   double eksum_new;
   double pnhf_pre;
   double scale;
   double delp;
   double randfor;

   dpress = (press - bath::atmsph) / units::prescon;
   delp = dpress;
   delpr = delp - 0.25 * (2.0 * eksum / vbox) / 3.0;
   pwinv = 1.0 / (3 * hmass_lp);
   hdot_old = hdot_lp; // h vel at -1/2

   pnhf_pre = pnhf_lp;
   pnhv_pre_lp = pnhv_lp;

   darray::copy(g::q0, n, vx, leapfrog_vxold);
   darray::copy(g::q0, n, vy, leapfrog_vyold);
   darray::copy(g::q0, n, vz, leapfrog_vzold);

   double xtlabc_old = xtlabc;
   hdot_old = hdot_lp;
   for (int i = 1; i <= 3; ++i) {
      xtlabc = xtlabc_old;
      hdot_lp = hdot_old;

      delpx = 3.0 * delp / xtlabc;
      h_random_force = normal<double>(0, prfwd);
      hdot_pre = hdot_lp;
      randfor = pbfact * h_random_force / dt;
      hdot_lp = palpha * hdot_lp + pwinv * delpx * vbox * pbfact / dt +
         randfor;                            // h vel at 1/2
      hdot_mid = 0.5 * (hdot_lp + hdot_pre); // h vel at 0

      sf = xtlabc + hdot_lp * dt;   // h at 1
      sh = sf - 0.5 * hdot_lp * dt; // h at 1/2
      // thermostat
      pnhf_lp = (2.0 * eksum_mid - df * kbt) / pnhm_lp;
      if (pnhf_pre == 0.0)
         pnhf_pre = pnhf_lp;
      pnhv_lp = pnhv_pre_lp + 0.5 * (pnhf_pre + pnhf_lp) * dt;

      eksum_new = 0.0;
      eksum = 0.0;
      eksum_mid = 0.0;
      scale = -(hdot_mid / xtlabc + pnhv_lp) * dt;
      // vx = vxnew + scale * (vold + vx)/2
      propagate_velocity_lp(vx, vy, vz, leapfrog_vxold, leapfrog_vyold,
                            leapfrog_vzold, leapfrog_vx, leapfrog_vy,
                            leapfrog_vz, scale, eksum_new,
                            eksum_mid); // eksum_new -> vx; eksum -> (vx+vold)/2


      // x = vx * dt + xold + scale * (x+xold)
      scale = 0.5 * hdot_lp / sh * dt;
      propagate_pos_lp(dt, xpos, ypos, zpos, vx, vy, vz, leapfrog_x, leapfrog_y,
                       leapfrog_z, scale);

      xtlabc = sf;
      delp = delpr + 0.25 * (2.0 * eksum_new / vbox) / 3.0;
   } // end of loop

   const bool userat = use_rattle();
   if (userat) {

      shake(dt, xpos, ypos, zpos, leapfrog_x, leapfrog_y, leapfrog_z);

      scale = 0.5 * hdot_lp * dt / sh;

      // xold = xold + scale * (x + xold)
      propagate_pos_lp2(dt, xpos, ypos, zpos, leapfrog_x, leapfrog_y,
                        leapfrog_z, scale);
      shake(dt, leapfrog_x, leapfrog_y, leapfrog_z, xpos, ypos, zpos);

      // vx =  (x - xold) / dt
      propagate_velocity_lp2(dt, vx, vy, vz, xpos, ypos, zpos, leapfrog_x,
                             leapfrog_y, leapfrog_z);


   } // userat

   darray::copy(g::q0, n, leapfrog_vxold, vx);
   darray::copy(g::q0, n, leapfrog_vyold, vy);
   darray::copy(g::q0, n, leapfrog_vzold, vz);
   // vx = 0.5 * (vxnew + vxold)
   eksum_new = 0.0;
   propagate_velocity_lp3(vx, vy, vz, leapfrog_vx, leapfrog_vy, leapfrog_vz,
                          leapfrog_vxold, leapfrog_vyold, leapfrog_vzold,
                          eksum_new); // eksum_new => vx

   pnhf_lp = (2.0 * eksum_new - df * kbt) / pnhm_lp;
   pnhv_lp = pnhv_pre_lp + 0.5 * (pnhf_pre + pnhf_lp) * dt;


   press = delp * units::prescon + bath::atmsph;
   lvec1.x = xtlabc;
   lvec2.y = xtlabc;
   lvec3.z = xtlabc;
   set_default_recip_box();

   pnh_lp = pnh_lp + pnhv_lp * dt + 0.5 * dt * dt * pnhf_lp;

   return;
}

void lf_lpiston_npt(int istep, time_prec dt_ps)
{
   int vers1 = rc_flag & calc::vmask;
   bool save = !(istep % inform::iwrite);
   double eksum_new = 0.0;

   if (!save)
      vers1 &= ~calc::energy;
   time_prec t2 = 0.5 * dt_ps;
   T_prec temp;
   const bool userat = use_rattle();
   bool com_done = (((istep - 1) % mdstuf::irest) == 0);
   if (istep == 1 || com_done) {

      if (istep == 1) {
         darray::copy(g::q0, n, leapfrog_x, xpos);
         darray::copy(g::q0, n, leapfrog_y, ypos);
         darray::copy(g::q0, n, leapfrog_z, zpos);

         if (userat) {
            shake(dt_ps, xpos, ypos, zpos, leapfrog_x, leapfrog_y, leapfrog_z);
            copy_pos_to_xyz(true);
            darray::copy(g::q0, n, leapfrog_x, xpos);
            darray::copy(g::q0, n, leapfrog_y, ypos);
            darray::copy(g::q0, n, leapfrog_z, zpos);
            energy(vers1);
         }

      } else {

         darray::copy(g::q0, n, xpos, leapfrog_x);
         darray::copy(g::q0, n, ypos, leapfrog_y);
         darray::copy(g::q0, n, zpos, leapfrog_z);
         if (userat) {
            shake(dt_ps, xpos, ypos, zpos, leapfrog_x, leapfrog_y, leapfrog_z);
            darray::copy(g::q0, n, leapfrog_x, xpos);
            darray::copy(g::q0, n, leapfrog_y, ypos);
            darray::copy(g::q0, n, leapfrog_z, zpos);
         }

         copy_pos_to_xyz(true);
         energy(vers1);
      }

      // in step 1, vx(new) and vxold are opposite from notation
      // v(-0.5) = v(0) - a(0)*dt/2
      propagate_velocity(-t2, leapfrog_vx, leapfrog_vy, leapfrog_vz, vx, vy, vz,
                         gx, gy, gz);
      // v(0.5) = v(0) + a(0)*dt/2
      propagate_velocity(t2, leapfrog_vxold, leapfrog_vyold, leapfrog_vzold, vx,
                         vy, vz, gx, gy, gz);

      if (userat) {

         // x = xold - vnew * dt
         propagate_pos_lf(-dt_ps, xpos, ypos, zpos, leapfrog_x, leapfrog_y,
                          leapfrog_z, leapfrog_vx, leapfrog_vy, leapfrog_vz);
         shake(dt_ps, xpos, ypos, zpos, leapfrog_x, leapfrog_y, leapfrog_z);
         // vnew =  (xold - x) / dt
         propagate_velocity_lp2(dt_ps, leapfrog_vx, leapfrog_vy, leapfrog_vz,
                                leapfrog_x, leapfrog_y, leapfrog_z, xpos, ypos,
                                zpos);

         // x = xold + vold * dt
         propagate_pos_lf(dt_ps, xpos, ypos, zpos, leapfrog_x, leapfrog_y,
                          leapfrog_z, leapfrog_vxold, leapfrog_vyold,
                          leapfrog_vzold);

         // propagate vx = vnew + a(1)*dt
         propagate_velocity(dt_ps, vx, vy, vz, leapfrog_vx, leapfrog_vy,
                            leapfrog_vz, gx, gy, gz);

         shake(dt_ps, xpos, ypos, zpos, leapfrog_x, leapfrog_y, leapfrog_z);

         // vold =  (x - xold) / dt
         propagate_velocity_lp2(dt_ps, leapfrog_vxold, leapfrog_vyold,
                                leapfrog_vzold, xpos, ypos, zpos, leapfrog_x,
                                leapfrog_y, leapfrog_z);
         // virshk
         shake2(dt_ps, vx, vy, vz, leapfrog_vxold, leapfrog_vyold,
                leapfrog_vzold, leapfrog_x, leapfrog_y, leapfrog_z);
         // vx = (vnew + vold) / 2
         eksum_new = 0.0;
         propagate_velocity_lp3(vx, vy, vz, leapfrog_vx, leapfrog_vy,
                                leapfrog_vz, leapfrog_vxold, leapfrog_vyold,
                                leapfrog_vzold, eksum_new); // eksum_new -> vx
      }

      // calculate eksum (at n+1/2) and eksum_old (at n-1/2)
      kinetic_leapfrog(temp);
      double vbox = volbox();
      double factor = units::prescon / vbox;
      double stress[3][3];

      for (int i = 0; i < 3; ++i) {
         for (int j = 0; j < 3; ++j) {
            stress[i][j] = factor * (-vir[3 * i + j]);
         }
      }
      press = (stress[0][0] + stress[1][1] + stress[2][2]) / 3;
      press += (eksum + eksum_old) * factor / 3;
   }

   // propagate x(0 -> 1): x = xold + v(0.5)*dt
   propagate_pos_lf(dt_ps, xpos, ypos, zpos, leapfrog_x, leapfrog_y, leapfrog_z,
                    leapfrog_vxold, leapfrog_vyold, leapfrog_vzold);
   lf_langevin_piston(dt_ps, press);
   copy_pos_to_xyz(true);
   energy(vers1);

   darray::copy(g::q0, n, leapfrog_x, xpos);
   darray::copy(g::q0, n, leapfrog_y, ypos);
   darray::copy(g::q0, n, leapfrog_z, zpos);

   // propagate v(0.5 -> 1.5 old): v = vold + a(1)*dt
   propagate_velocity(dt_ps, leapfrog_vx, leapfrog_vy, leapfrog_vz,
                      leapfrog_vxold, leapfrog_vyold, leapfrog_vzold, gx, gy,
                      gz);

   if (userat) {
      darray::copy(g::q0, n, vx, leapfrog_vx);
      darray::copy(g::q0, n, vy, leapfrog_vy);
      darray::copy(g::q0, n, vz, leapfrog_vz);

      // x = xold + vnew * dt
      propagate_pos_lf(dt_ps, xpos, ypos, zpos, leapfrog_x, leapfrog_y,
                       leapfrog_z, leapfrog_vx, leapfrog_vy, leapfrog_vz);

      shake(dt_ps, xpos, ypos, zpos, leapfrog_x, leapfrog_y, leapfrog_z);

      // vnew =  (x - xold) / dt
      propagate_velocity_lp2(dt_ps, leapfrog_vx, leapfrog_vy, leapfrog_vz, xpos,
                             ypos, zpos, leapfrog_x, leapfrog_y, leapfrog_z);

      // do virial correction
      shake2(dt_ps, vx, vy, vz, leapfrog_vx, leapfrog_vy, leapfrog_vz,
             leapfrog_x, leapfrog_y, leapfrog_z);
   }

   // vx = (vnew + vold) / 2
   eksum_new = 0.0;
   eksum_mid = 0.0;
   propagate_velocity_lp3(vx, vy, vz, leapfrog_vx, leapfrog_vy, leapfrog_vz,
                          leapfrog_vxold, leapfrog_vyold, leapfrog_vzold,
                          eksum_mid);


   kinetic_leapfrog(temp); // get eksum and eksum_old

   double vbox = volbox();
   double factor = units::prescon / vbox;
   double stress[3][3];
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         stress[i][j] = factor * (-vir[3 * i + j]);
      }
   }
   press = (stress[0][0] + stress[1][1] + stress[2][2]) / 3;
   press += (eksum + eksum_old) * factor / 3;
   swap_velocity(leapfrog_vx, leapfrog_vy, leapfrog_vz, leapfrog_vxold,
                 leapfrog_vyold, leapfrog_vzold);
}


void kinetic_leapfrog(T_prec& temp)
{
   // Ek at +1/2
   T_prec t1;
   energy_prec ekin1[3][3];
   kinetic_explicit(t1, eksum, ekin1, leapfrog_vx, leapfrog_vy, leapfrog_vz);


   // Ek at -1/2
   T_prec t2;
   energy_prec ekin2[3][3];
   kinetic_explicit(t2, eksum_old, ekin2, leapfrog_vxold, leapfrog_vyold,
                    leapfrog_vzold);


   ekin[0][0] = 0.5 * (ekin1[0][0] + ekin2[0][0]);
   ekin[0][1] = 0.5 * (ekin1[0][1] + ekin2[0][1]);
   ekin[0][2] = 0.5 * (ekin1[0][2] + ekin2[0][2]);
   ekin[1][0] = 0.5 * (ekin1[1][0] + ekin2[1][0]);
   ekin[1][1] = 0.5 * (ekin1[1][1] + ekin2[1][1]);
   ekin[1][2] = 0.5 * (ekin1[1][2] + ekin2[1][2]);
   ekin[2][0] = 0.5 * (ekin1[2][0] + ekin2[2][0]);
   ekin[2][1] = 0.5 * (ekin1[2][1] + ekin2[2][1]);
   ekin[2][2] = 0.5 * (ekin1[2][2] + ekin2[2][2]);
   temp = 0.5 * (t1 + t2);
}


void shake2(time_prec dt, const vel_prec* vxold, const vel_prec* vyold,
            const vel_prec* vzold, const vel_prec* vxnew, const vel_prec* vynew,
            const vel_prec* vznew, const pos_prec* xold, const pos_prec* yold,
            const pos_prec* zold)
{
   darray::zero(g::q0, buffer_size(), vir_buf);


   shake2_acc(dt, vxold, vyold, vzold, vxnew, vynew, vznew, xold, yold, zold);


   virial_prec v[9];
   virial_reduce(v, vir_buf);
   for (int iv = 0; iv < 9; ++iv) {
      vir[iv] += v[iv];
   }
}


//====================================================================//


void swap_velocity(vel_prec* vxnew, vel_prec* vynew, vel_prec* vznew,
                   vel_prec* vxold, vel_prec* vyold, vel_prec* vzold)
{
   swap_velocity_acc(vxnew, vynew, vznew, vxold, vyold, vzold);
}


void propagate_pos_lp(time_prec dt, pos_prec* x_lp, pos_prec* y_lp,
                      pos_prec* z_lp, const vel_prec* vx_lp,
                      const vel_prec* vy_lp, const vel_prec* vz_lp,
                      const pos_prec* xold_lp, const pos_prec* yold_lp,
                      const pos_prec* zold_lp, double scale)
{

   propagate_pos_lp_acc(dt, x_lp, y_lp, z_lp, vx_lp, vy_lp, vz_lp, xold_lp,
                        yold_lp, zold_lp, scale);
}


void propagate_pos_lp2(time_prec dt, const pos_prec* x_lp, const pos_prec* y_lp,
                       const pos_prec* z_lp, pos_prec* xold_lp,
                       pos_prec* yold_lp, pos_prec* zold_lp, double scale)
{
   propagate_pos_lp2_acc(dt, x_lp, y_lp, z_lp, xold_lp, yold_lp, zold_lp,
                         scale);
}


void propagate_pos_lf(time_prec dt, pos_prec* qx, pos_prec* qy, pos_prec* qz,
                      const pos_prec* qxold, const pos_prec* qyold,
                      const pos_prec* qzold, const vel_prec* vlx,
                      const vel_prec* vly, const vel_prec* vlz)
{
   propagate_pos_lf_acc(dt, qx, qy, qz, qxold, qyold, qzold, vlx, vly, vlz);
}


void propagate_velocity_lp(vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp,
                           const vel_prec* vxnew_lp, const vel_prec* vynew_lp,
                           const vel_prec* vznew_lp, const vel_prec* vxold_lp,
                           const vel_prec* vyold_lp, const vel_prec* vzold_lp,
                           const double scale, energy_prec& eksum_new,
                           energy_prec& eksum_old)
{
   propagate_velocity_lp_acc(vx_lp, vy_lp, vz_lp, vxnew_lp, vynew_lp, vznew_lp,
                             vxold_lp, vyold_lp, vzold_lp, scale, eksum_new,
                             eksum_old);
}


void propagate_velocity_lp2(time_prec dt, vel_prec* vx_lp, vel_prec* vy_lp,
                            vel_prec* vz_lp, const pos_prec* x_lp,
                            const pos_prec* y_lp, const pos_prec* z_lp,
                            const pos_prec* xold_lp, const pos_prec* yold_lp,
                            const pos_prec* zold_lp)
{
   propagate_velocity_lp2_acc(dt, vx_lp, vy_lp, vz_lp, x_lp, y_lp, z_lp,
                              xold_lp, yold_lp, zold_lp);
}


void propagate_velocity_lp3(vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp,
                            const vel_prec* vxnew_lp, const vel_prec* vynew_lp,
                            const vel_prec* vznew_lp, const vel_prec* vxold_lp,
                            const vel_prec* vyold_lp, const vel_prec* vzold_lp,
                            energy_prec& eksum_new)
{
   propagate_velocity_lp3_acc(vx_lp, vy_lp, vz_lp, vxnew_lp, vynew_lp, vznew_lp,
                              vxold_lp, vyold_lp, vzold_lp, eksum_new);
}
}
