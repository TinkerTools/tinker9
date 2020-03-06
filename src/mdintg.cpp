#define TINKER_ENABLE_LOG 0
#include "log.h"


#include "darray.h"
#include "energy.h"
#include "io_fort_str.h"
#include "md_calc.h"
#include "md_egv.h"
#include "md_pt.h"
#include "md_save.h"
#include "mdintg.h"
#include "mdpq.h"
#include <cassert>
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>


TINKER_NAMESPACE_BEGIN
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
      TINKER_LOG("Integrating Step {:10d}", istep);
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
      if (intg == respa_fast_slow)
         darray::deallocate(gx1, gy1, gz1, gx2, gy2, gz2);

      if (barostat == MONTE_CARLO_BAROSTAT)
         darray::deallocate(x_pmonte, y_pmonte, z_pmonte);

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
         else if (br == "MONTECARLO") {
            barostat = MONTE_CARLO_BAROSTAT;
            darray::allocate(n, &x_pmonte, &y_pmonte, &z_pmonte);
         } else
            assert(false);
      } else {
         barostat = NONE_BAROSTAT;
      }

      // Only gradient is necessary to start a simulation.
      fstr_view itg = mdstuf::integrate;
      intg = nullptr;
      if (itg == "VERLET") {
         intg = velocity_verlet;
         // need full gradient to start/restart the simulation
         energy(calc::grad);
      } else if (itg == "STOCHASTIC") {
      } else if (itg == "BAOAB") {
      } else if (itg == "BUSSI") {
      } else if (itg == "NOSE-HOOVER") {
      } else if (itg == "GHMC") {
      } else if (itg == "RIGIDBODY") {
      } else if (itg == "RESPA") {
         intg = respa_fast_slow;
         // need fast and slow gradients to start/restart the simulation
         darray::allocate(n, &gx1, &gy1, &gz1, &gx2, &gy2, &gz2);
         // save fast gradients to gx1 etc.
         energy(calc::grad, RESPA_FAST, respa_tsconfig());
         darray::copy(PROCEED_NEW_Q, n, gx1, gx);
         darray::copy(PROCEED_NEW_Q, n, gy1, gy);
         darray::copy(PROCEED_NEW_Q, n, gz1, gz);

         // save slow gradients to gx2 etc.
         energy(calc::grad, RESPA_SLOW, respa_tsconfig());
         darray::copy(PROCEED_NEW_Q, n, gx2, gx);
         darray::copy(PROCEED_NEW_Q, n, gy2, gy);
         darray::copy(PROCEED_NEW_Q, n, gz2, gz);
      } else {
         // beeman
      }
      assert(intg);
   }
}
TINKER_NAMESPACE_END
