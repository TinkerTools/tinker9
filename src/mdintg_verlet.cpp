#include "box.h"
#include "energy.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include "mdpt.h"
#include "rattle.h"
#include <chrono>
#include <random>
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/stodyn.hh>
#include <tinker/detail/units.hh>


namespace tinker {
void velocity_verlet(int istep, time_prec dt_ps)
{
   int vers0 = rc_flag & calc::vmask;
   int vers1 = vers0;

   bool save = !(istep % inform::iwrite);
   bool mcbaro = do_pmonte;
   if (barostat == MONTE_CARLO_BAROSTAT) {
      // toggle off the calc::virial bit if Monte Carlo Barostat is in use
      vers1 &= ~calc::virial;
   }
   // toggle off the calc::energy bit if neither save nor mcbaro
   if (!save && !mcbaro)
      vers1 &= ~calc::energy;

   time_prec dt_2 = 0.5f * dt_ps;

   // gradients were calculated in integrate_data()
   // v += a * dt/2
   propagate_velocity(dt_2, gx, gy, gz);

   const bool userat = use_rattle();
   if (userat) {
      darray::copy(PROCEED_NEW_Q, n, rattle_xold, xpos);
      darray::copy(PROCEED_NEW_Q, n, rattle_yold, ypos);
      darray::copy(PROCEED_NEW_Q, n, rattle_zold, zpos);
   }
   // s += v * dt
   propagate_pos(dt_ps);
   if (userat)
      rattle(dt_ps, rattle_xold, rattle_yold, rattle_zold);
   copy_pos_to_xyz(true);

   // update gradient
   energy(vers1);

   // half-step corrections for certain thermostats and barostats
   T_prec temp;
   temper2(dt_ps, temp);
   pressure2(esum, temp);

   // v += a * dt/2
   propagate_velocity(dt_2, gx, gy, gz);
   if (userat)
      rattle2(dt_ps, vers1 bitand calc::virial);

   // full-step corrections
   temper(dt_ps, temp);
   pressure();
}


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


void langevin_piston(time_prec dt, virial_prec press)
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
   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   std::default_random_engine generator(seed);
   std::normal_distribution<double> random_gaussian(0.0, prfwd);
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

   darray::copy(PROCEED_NEW_Q, n, vx, leapfrog_vxold);
   darray::copy(PROCEED_NEW_Q, n, vy, leapfrog_vyold);
   darray::copy(PROCEED_NEW_Q, n, vz, leapfrog_vzold);

   double xtlabc_old = xtlabc;
   hdot_old = hdot_lp;
   for (int i = 1; i <= 3; ++i) {
      xtlabc = xtlabc_old;
      hdot_lp = hdot_old;

      delpx = 3.0 * delp / xtlabc;
      h_random_force = random_gaussian(generator);
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

   darray::copy(PROCEED_NEW_Q, n, leapfrog_vxold, vx);
   darray::copy(PROCEED_NEW_Q, n, leapfrog_vyold, vy);
   darray::copy(PROCEED_NEW_Q, n, leapfrog_vzold, vz);
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

void lpiston_npt(int istep, time_prec dt_ps)
{
   int vers1 = rc_flag & calc::vmask;
   bool save = !(!istep % inform::iwrite);
   double eksum_new = 0.0;

   if (!save)
      vers1 &= ~calc::energy;
   time_prec t2 = 0.5 * dt_ps;
   T_prec temp;
   const bool userat = use_rattle();
   bool com_done = (((istep - 1) % mdstuf::irest) == 0);
   if (istep == 1 || com_done) {

      if (istep == 1) {
         darray::copy(PROCEED_NEW_Q, n, leapfrog_x, xpos);
         darray::copy(PROCEED_NEW_Q, n, leapfrog_y, ypos);
         darray::copy(PROCEED_NEW_Q, n, leapfrog_z, zpos);

         if (userat) {
            shake(dt_ps, xpos, ypos, zpos, leapfrog_x, leapfrog_y, leapfrog_z);
            copy_pos_to_xyz(true);
            darray::copy(PROCEED_NEW_Q, n, leapfrog_x, xpos);
            darray::copy(PROCEED_NEW_Q, n, leapfrog_y, ypos);
            darray::copy(PROCEED_NEW_Q, n, leapfrog_z, zpos);
            energy(vers1);
         }

      } else {

         darray::copy(PROCEED_NEW_Q, n, xpos, leapfrog_x);
         darray::copy(PROCEED_NEW_Q, n, ypos, leapfrog_y);
         darray::copy(PROCEED_NEW_Q, n, zpos, leapfrog_z);
         if (userat) {
            shake(dt_ps, xpos, ypos, zpos, leapfrog_x, leapfrog_y, leapfrog_z);
            darray::copy(PROCEED_NEW_Q, n, leapfrog_x, xpos);
            darray::copy(PROCEED_NEW_Q, n, leapfrog_y, ypos);
            darray::copy(PROCEED_NEW_Q, n, leapfrog_z, zpos);
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
      temper_leapfrog(dt_ps, temp);
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
   langevin_piston(dt_ps, press);
   copy_pos_to_xyz(true);
   energy(vers1);

   darray::copy(PROCEED_NEW_Q, n, leapfrog_x, xpos);
   darray::copy(PROCEED_NEW_Q, n, leapfrog_y, ypos);
   darray::copy(PROCEED_NEW_Q, n, leapfrog_z, zpos);

   // propagate v(0.5 -> 1.5 old): v = vold + a(1)*dt
   propagate_velocity(dt_ps, leapfrog_vx, leapfrog_vy, leapfrog_vz,
                      leapfrog_vxold, leapfrog_vyold, leapfrog_vzold, gx, gy,
                      gz);

   if (userat) {
      darray::copy(PROCEED_NEW_Q, n, vx, leapfrog_vx);
      darray::copy(PROCEED_NEW_Q, n, vy, leapfrog_vy);
      darray::copy(PROCEED_NEW_Q, n, vz, leapfrog_vz);

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


   temper_leapfrog(dt_ps, temp); // get eksum and eksum_old

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
}
