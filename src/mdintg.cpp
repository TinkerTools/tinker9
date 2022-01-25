#define TINKER_ENABLE_LOG 0
#include "tool/log.h"

#include "energy.h"
#include "intg/enum.h"
#include "intg/intgBasic.h"
#include "intg/intgLeapFrogLP.h"
#include "intg/intgNhc96.h"
#include "intg/intgRespa.h"
#include "intg/intgVerlet.h"
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
void mdrest_acc(int istep);
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

static BasicIntegrator* intg;

void propagate(int nsteps, time_prec dt_ps)
{
   for (int istep = 1; istep <= nsteps; ++istep) {
      TINKER_LOG("Integrating Step %10d", istep);
      intg->dynamic(istep, dt_ps);

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
      delete intg;
      intg = nullptr;
   }

   if (op & rc_init) {
      ThermostatEnum thermostat = ThermostatEnum::Null;
      if (bath::isothermal) {
         fstr_view th = bath::thermostat;
         if (th == "ANDERSEN")
            thermostat = ThermostatEnum::Andersen;
         else if (th == "BERENDSEN")
            thermostat = ThermostatEnum::Berendsen;
         else if (th == "BUSSI")
            thermostat = ThermostatEnum::Bussi;
         else if (th == "LPISTON")
            thermostat = ThermostatEnum::LeapFrogLP;
         else if (th == "NOSE-HOOVER")
            thermostat = ThermostatEnum::Nhc96;
         else
            assert(false);
      }

      BarostatEnum barostat = BarostatEnum::Null;
      if (bath::isobaric) {
         fstr_view br = bath::barostat;
         if (br == "BERENDSEN")
            barostat = BarostatEnum::Berendsen;
         else if (br == "BUSSI")
            barostat = BarostatEnum::Bussi;
         else if (br == "LANGEVIN")
            barostat = BarostatEnum::Langevin;
         else if (br == "LPISTON")
            barostat = BarostatEnum::LeapFrogLP;
         else if (br == "MONTECARLO")
            barostat = BarostatEnum::MonteCarlo;
         else if (br == "NOSE-HOOVER")
            barostat = BarostatEnum::Nhc96;
         else
            assert(false);
      }

      IntegratorEnum integrator = IntegratorEnum::Beeman;
      fstr_view itg = mdstuf::integrate;
      if (itg == "VERLET") {
         integrator = IntegratorEnum::Verlet;
      } else if (itg == "LPISTON") {
         integrator = IntegratorEnum::LeapFrogLP;
         thermostat = ThermostatEnum::LeapFrogLP;
         barostat = BarostatEnum::LeapFrogLP;
      } else if (itg == "NOSE-HOOVER") {
         integrator = IntegratorEnum::Nhc96;
         thermostat = ThermostatEnum::Nhc96;
         barostat = BarostatEnum::Nhc96;
      } else if (itg == "RESPA") {
         integrator = IntegratorEnum::Respa;
      }

      if (thermostat == ThermostatEnum::LeapFrogLP and
          barostat == BarostatEnum::LeapFrogLP)
         integrator = IntegratorEnum::LeapFrogLP;
      else if (barostat == BarostatEnum::Langevin) {
         if (itg == "VERLET" or itg == "RESPA")
            integrator = IntegratorEnum::LangevinNpt;
      } else if (thermostat == ThermostatEnum::Nhc96 and
                 barostat == BarostatEnum::Nhc96)
         integrator = IntegratorEnum::Nhc96;

      intg = nullptr;
      if (integrator == IntegratorEnum::LangevinNpt) {
      } else if (integrator == IntegratorEnum::LeapFrogLP)
         intg = new LeapFrogLPIntegrator;
      else if (integrator == IntegratorEnum::Nhc96)
         intg = new Nhc96Integrator;
      else if (integrator == IntegratorEnum::Respa)
         intg = new RespaIntegrator(thermostat, barostat);
      else if (integrator == IntegratorEnum::Verlet)
         intg = new VerletIntegrator(thermostat, barostat);
      else if (integrator == IntegratorEnum::Beeman)
         TINKER_THROW("Beeman integrator is not available.");
   }
}
}
