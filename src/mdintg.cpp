#define TINKER_ENABLE_LOG 0
#include "tool/log.h"


#include "energy.h"
#include "lf_lpiston.h"
#include "lpiston.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include "mdpt.h"
#include "mdsave.h"
#include "nose.h"
#include "random.h"
#include "rattle.h"
#include "tool/darray.h"
#include "tool/io_fort_str.h"
#include <cassert>
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/stodyn.hh>
#include <tinker/detail/units.hh>


namespace tinker {
void mdrest(int istep)
{
   mdrest_acc(istep);
}


void md_data(rc_op op)
{
   if ((calc::md & rc_flag) == 0)
      return;

   rc_man intg42_{integrate_data, op};
   rc_man save42_{mdsave_data, op};
}


//====================================================================//


namespace {
void (*intg)(int, time_prec);
}


void propagate(int nsteps, time_prec dt_ps)
{
   for (int istep = 1; istep <= nsteps; ++istep) {
      do_pmonte = false;
      if (barostat == MONTE_CARLO_BAROSTAT) {
         double rdm = random<double>();
         if (rdm < 1.0 / bath::voltrial)
            do_pmonte = true;
      }

      TINKER_LOG("Integrating Step %10d", istep);
      intg(istep, dt_ps);

      // mdstat
      bool save = (istep % inform::iwrite == 0);
      if (save || (istep % BOUNDS_EVERY_X_STEPS) == 0)
         bounds();
      if (save) {
         T_prec temp;
         kinetic(temp);
         mdsave_async(istep, dt_ps);
      }
      mdrest(istep);
   }
   mdsave_synchronize();
}


void integrate_data(rc_op op)
{
   if (op & rc_dealloc) {
      if (intg == lf_lpiston_npt) {
         darray::deallocate(leapfrog_x, leapfrog_y, leapfrog_z);
         darray::deallocate(leapfrog_vx, leapfrog_vy, leapfrog_vz,
                            leapfrog_vxold, leapfrog_vyold, leapfrog_vzold);
      }

      if (intg == vv_lpiston_npt) {
         if (use_rattle())
            darray::deallocate(lp_molpres_buf);
      }

      if (intg == respa_fast_slow)
         darray::deallocate(gx1, gy1, gz1, gx2, gy2, gz2);

      if (barostat == MONTE_CARLO_BAROSTAT) {
         darray::deallocate(x_pmonte, y_pmonte, z_pmonte);
         darray::deallocate(vx_pmonte, vy_pmonte, vz_pmonte);
      }

      intg = nullptr;
   }

   if (op & rc_init) {
      if (bath::isothermal) {
         fstr_view th = bath::thermostat;
         if (th == "BERENDSEN")
            thermostat = BERENDSEN_THERMOSTAT;
         else if (th == "BUSSI")
            thermostat = BUSSI_THERMOSTAT;
         else if (th == "ANDERSEN")
            thermostat = ANDERSEN_THERMOSTAT;
         else if (th == "NOSE-HOOVER")
            thermostat = NOSE_HOOVER_CHAIN_THERMOSTAT;
         else if (th == "LPISTON")
            thermostat = LANGEVIN_PISTON_THERMOSTAT;
         else if (th == "VVLP")
            thermostat = VV_LPISTON_THERMOSTAT;
         else
            assert(false);
      } else {
         thermostat = NONE_THERMOSTAT;
      }

      if (bath::isobaric) {
         fstr_view br = bath::barostat;
         if (br == "BERENDSEN")
            barostat = BERENDSEN_BAROSTAT;
         else if (br == "BUSSI")
            barostat = BUSSI_BAROSTAT;
         else if (br == "NOSE-HOOVER")
            barostat = NOSE_HOOVER_CHAIN_BAROSTAT;
         else if (br == "LPISTON")
            barostat = LANGEVIN_PISTON_BAROSTAT;
         else if (br == "VVLP")
            barostat = VV_LPISTON_BAROSTAT;
         else if (br == "MONTECARLO") {
            barostat = MONTE_CARLO_BAROSTAT;
            darray::allocate(n, &x_pmonte, &y_pmonte, &z_pmonte);
            darray::allocate(n, &vx_pmonte, &vy_pmonte, &vz_pmonte);
         } else
            assert(false);
      } else {
         barostat = NONE_BAROSTAT;
      }

      fstr_view itg = mdstuf::integrate;
      intg = nullptr;
      if (itg == "VERLET") {
         intg = velocity_verlet;
      } else if (itg == "LPISTON") {
         intg = lf_lpiston_npt;
         thermostat = LANGEVIN_PISTON_THERMOSTAT;
         barostat = LANGEVIN_PISTON_BAROSTAT;
      } else if (itg == "VVLP") {
         intg = vv_lpiston_npt;
         thermostat = VV_LPISTON_THERMOSTAT;
         barostat = VV_LPISTON_BAROSTAT;
      } else if (itg == "NOSE-HOOVER") {
         intg = nhc_npt;
         thermostat = NOSE_HOOVER_CHAIN_THERMOSTAT;
         barostat = NOSE_HOOVER_CHAIN_BAROSTAT;
      } else if (itg == "RESPA") {
         intg = respa_fast_slow;
      }

      if (thermostat == LANGEVIN_PISTON_THERMOSTAT and
          barostat == LANGEVIN_PISTON_BAROSTAT) {
         intg = lf_lpiston_npt;
      } else if (thermostat == VV_LPISTON_THERMOSTAT and
                 barostat == VV_LPISTON_BAROSTAT) {
         intg = vv_lpiston_npt;
      } else if (thermostat == NOSE_HOOVER_CHAIN_THERMOSTAT and
                 barostat == NOSE_HOOVER_CHAIN_BAROSTAT) {
         intg = nhc_npt;
      }

      // Only gradient is necessary to start a simulation.
      if (intg == velocity_verlet) {
         // need full gradient to start/restart the simulation
         energy(calc::grad);
      } else if (intg == lf_lpiston_npt) {
         darray::allocate(n, &leapfrog_x, &leapfrog_y, &leapfrog_z);
         darray::allocate(n, &leapfrog_vx, &leapfrog_vy, &leapfrog_vz,
                          &leapfrog_vxold, &leapfrog_vyold, &leapfrog_vzold);
         energy(calc::v1);
      } else if (intg == vv_lpiston_npt) {
         double ekt = units::gasconst * bath::kelvin;
         vbar = 0;
         gbar = 0;
         qbar = (mdstuf::nfree + 3) * ekt * bath::taupres * bath::taupres;
         for (int i = 0; i < maxnose; ++i) {
            vnh[i] = 0;
            gnh[i] = 0;
            qnh[i] = ekt * bath::tautemp * bath::tautemp;
         }
         qnh[0] = mdstuf::nfree * ekt * bath::tautemp * bath::tautemp;
         lp_alpha = 1.0;
         if (n > 1)
            lp_alpha = 1.0 + 1.0 / (n - 1);
         energy(calc::v6);
         if (use_rattle()) {
            darray::allocate(buffer_size(), &lp_molpres_buf);
            lp_molpressure(lp_alpha, lp_molpres);
         }

         printf("\n");
         printf(" Friction                        %12.4lf /ps\n",
                stodyn::friction);
         printf(" Time constant for the const-T   %12.4lf ps\n", bath::tautemp);
         printf(" Time constant for the const-P   %12.4lf ps\n", bath::taupres);
         printf("\n");
      } else if (intg == nhc_npt) {
         if (use_rattle()) {
            TINKER_THROW(
               "Constraints under NH-NPT require the ROLL algorithm.");
         }
         double ekt = units::gasconst * bath::kelvin;
         vbar = 0;
         qbar = (mdstuf::nfree + 1) * ekt * bath::taupres * bath::taupres;
         gbar = 0;
         for (int i = 0; i < maxnose; ++i) {
            vnh[i] = 0;
            qnh[i] = ekt * bath::tautemp * bath::tautemp;
            gnh[i] = 0;
         }
         qnh[0] *= mdstuf::nfree;
         energy(calc::v6);
      } else if (intg == respa_fast_slow) {
         // need fast and slow gradients to start/restart the simulation
         darray::allocate(n, &gx1, &gy1, &gz1, &gx2, &gy2, &gz2);

         // save fast gradients to gx1 etc.
         energy(calc::grad, RESPA_FAST, respa_tsconfig());
         darray::copy(g::q0, n, gx1, gx);
         darray::copy(g::q0, n, gy1, gy);
         darray::copy(g::q0, n, gz1, gz);

         // save slow gradients to gx2 etc.
         energy(calc::grad, RESPA_SLOW, respa_tsconfig());
         darray::copy(g::q0, n, gx2, gx);
         darray::copy(g::q0, n, gy2, gy);
         darray::copy(g::q0, n, gz2, gz);
      } else if (intg == nullptr) {
         // beeman
         TINKER_THROW("Beeman integrator is not available.");
      }
   }
}
}
