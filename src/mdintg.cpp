#define TINKER_ENABLE_LOG 0
#include "tool/log.h"


#include "energy.h"
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

   rc_man intg42{integrate_data, op};
   rc_man save42{mdsave_data, op};
}


//====================================================================//


namespace {
void (*intg)(int, time_prec);
}


void propagate(int nsteps, time_prec dt_ps)
{
   for (int istep = 1; istep <= nsteps; ++istep) {
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
      if (intg == vv_lpiston_npt) {
         vv_lpiston_destory();
      }

      

      if (barostat == MONTE_CARLO_BAROSTAT) {
         
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
            thermostat = LEAPFROG_LPISTON_THERMOSTAT;
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
            barostat = LEAPFROG_LPISTON_BAROSTAT;
         else if (br == "LANGEVIN")
            barostat = LANGEVIN_BAROSTAT;
         else if (br == "MONTECARLO") {
            barostat = MONTE_CARLO_BAROSTAT;
            
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
      } else if (itg == "NOSE-HOOVER") {
      } else if (itg == "RESPA") {
         intg = respa_fast_slow;
      }

      if (thermostat == LEAPFROG_LPISTON_THERMOSTAT and
          barostat == LEAPFROG_LPISTON_BAROSTAT) {
      } else if (barostat == LANGEVIN_BAROSTAT) {
         if (itg == "VERLET" or itg == "RESPA")
            intg = vv_lpiston_npt;
      } else if (thermostat == NOSE_HOOVER_CHAIN_THERMOSTAT and
                 barostat == NOSE_HOOVER_CHAIN_BAROSTAT) {
      }

      // Only gradient is necessary to start a simulation.
      if (intg == velocity_verlet) {
         // need full gradient to start/restart the simulation
         energy(calc::grad);
      } // else if (intg == lf_lpiston_npt) {
        // TODO
      // }
      else if (intg == vv_lpiston_npt) {
         vv_lpiston_init();
      } else if (intg == respa_fast_slow) {
      } else if (intg == nullptr) {
         // beeman
         TINKER_THROW("Beeman integrator is not available.");
      }
   }
}
}
