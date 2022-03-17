#define TINKER_ENABLE_LOG 0
#include "tool/log.h"

#include "energy.h"
#include "itgEnum.h"
#include "itgiBasic.h"
#include "itgiLP22.h"
#include "itgiLeapFrogLP.h"
#include "itgiNhc06.h"
#include "itgiNhc96.h"
#include "itgiRespa.h"
#include "itgiVerlet.h"
#include "lpiston.h"
#include "mathfunc_pow2.h"
#include "md.h"
#include "nose.h"
#include "random.h"
#include "rattle.h"
#include "tool/darray.h"
#include "tool/io.h"
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

void mdData(RcOp op)
{
   if ((calc::md & rc_flag) == 0)
      return;

   RcMan intg42{mdIntegrateData, op};
   RcMan save42{mdsaveData, op};
}

//====================================================================//

static BasicIntegrator* intg;

void mdPropagate(int nsteps, time_prec dt_ps)
{
   for (int istep = 1; istep <= nsteps; ++istep) {
      TINKER_LOG("Integrating Step %10d", istep);
      intg->dynamic(istep, dt_ps);

      // mdstat
      bool save = (istep % inform::iwrite == 0);
      if (save || (istep % BOUNDS_EVERY_X_STEPS) == 0)
         mdBounds();
      if (save) {
         T_prec temp;
         mdKinetic(temp);
         mdsaveAsync(istep, dt_ps);
      }
      mdrest(istep);
   }
   mdsaveSynchronize();
}

void mdIntegrateData(RcOp op)
{
   if (op & rc_dealloc) {
      delete intg;
      intg = nullptr;
   }

   if (op & rc_init) {
      ThermostatEnum thermostat = ThermostatEnum::Null;
      if (bath::isothermal) {
         FstrView th = bath::thermostat;
         if (th == "ANDERSEN")
            thermostat = ThermostatEnum::Andersen;
         else if (th == "BERENDSEN")
            thermostat = ThermostatEnum::Berendsen;
         else if (th == "BUSSI")
            thermostat = ThermostatEnum::Bussi;
         else if (th == "NOSE-HOOVER")
            thermostat = ThermostatEnum::Nhc;
         else if (th == "LPISTON")
            thermostat = ThermostatEnum::m_LeapFrogLP;
         else
            assert(false);
      }

      BarostatEnum barostat = BarostatEnum::Null;
      if (bath::isobaric) {
         FstrView br = bath::barostat;
         if (br == "BERENDSEN")
            barostat = BarostatEnum::Berendsen;
         else if (br == "BUSSI")
            barostat = BarostatEnum::Bussi;
         else if (br == "LANGEVIN")
            barostat = BarostatEnum::LP2022;
         else if (br == "MONTECARLO")
            barostat = BarostatEnum::MonteCarlo;
         else if (br == "NHC06")
            barostat = BarostatEnum::Nhc2006;
         else if (br == "LPISTON")
            barostat = BarostatEnum::m_LeapFrogLP;
         else if (br == "NOSE-HOOVER")
            barostat = BarostatEnum::m_Nhc1996;
         else
            assert(false);
      }

      IntegratorEnum integrator = IntegratorEnum::Beeman;
      if (thermostat == ThermostatEnum::m_LeapFrogLP and barostat == BarostatEnum::m_LeapFrogLP)
         integrator = IntegratorEnum::LeapFrogLP;
      else if (thermostat == ThermostatEnum::Nhc and barostat == BarostatEnum::m_Nhc1996)
         integrator = IntegratorEnum::Nhc1996;

      FstrView itg = mdstuf::integrate;
      if (itg == "RESPA") {
         integrator = IntegratorEnum::Respa;
      } else if (itg == "VERLET") {
         integrator = IntegratorEnum::Verlet;
      } else if (itg == "LPISTON") {
         integrator = IntegratorEnum::LeapFrogLP;
         thermostat = ThermostatEnum::m_LeapFrogLP;
         barostat = BarostatEnum::m_LeapFrogLP;
      } else if (itg == "NOSE-HOOVER") {
         integrator = IntegratorEnum::Nhc1996;
         thermostat = ThermostatEnum::m_Nhc1996;
         barostat = BarostatEnum::m_Nhc1996;
      }

      bool isNRespa1 = true;
      if (integrator == IntegratorEnum::Verlet or integrator == IntegratorEnum::Respa) {
         if (barostat == BarostatEnum::LP2022) {
            if (integrator == IntegratorEnum::Respa)
               isNRespa1 = false;
            integrator = IntegratorEnum::LP2022;
            thermostat = ThermostatEnum::m_LP2022;
         } else if (barostat == BarostatEnum::Nhc2006) {
            if (integrator == IntegratorEnum::Respa)
               isNRespa1 = false;
            integrator = IntegratorEnum::Nhc2006;
            thermostat = ThermostatEnum::m_Nhc2006;
         }
      }

      intg = nullptr;
      if (integrator == IntegratorEnum::Respa)
         intg = new RespaIntegrator(thermostat, barostat);
      else if (integrator == IntegratorEnum::Verlet)
         intg = new VerletIntegrator(thermostat, barostat);
      else if (integrator == IntegratorEnum::LeapFrogLP)
         intg = new LeapFrogLPIntegrator;
      else if (integrator == IntegratorEnum::LP2022)
         intg = new LP22Integrator(isNRespa1);
      else if (integrator == IntegratorEnum::Nhc1996)
         intg = new Nhc96Integrator;
      else if (integrator == IntegratorEnum::Nhc2006)
         intg = new Nhc06Integrator(isNRespa1);
      else if (integrator == IntegratorEnum::Beeman)
         TINKER_THROW("Beeman integrator is not available.");
      intg->printDetail(stdout);
   }
}

const TimeScaleConfig& mdRespaTsconfig()
{
   constexpr int fast = floor_log2_constexpr(RESPA_FAST); // short-range
   constexpr int slow = floor_log2_constexpr(RESPA_SLOW); // long-range
   static TimeScaleConfig tsconfig{
      {"ebond", fast},
      {"eangle", fast},
      {"estrbnd", fast},
      {"eurey", fast},
      {"eopbend", fast},
      {"etors", fast},
      {"eimprop", fast},
      {"eimptor", fast},
      {"epitors", fast},
      {"estrtor", fast},
      {"eangtor", fast},
      {"etortor", fast},
      {"egeom", fast},

      {"evalence", fast},

      {"evdw", slow},

      {"echarge", slow},
      {"echglj", slow},

      {"emplar", slow},
      {"empole", slow},
      {"epolar", slow},

      {"empole_chgpen", slow},
      {"epolar_chgpen", slow},

      {"echgtrn", slow},
      {"edisp", slow},
      {"erepel", slow},
      {"ehippo", slow},
   };
   return tsconfig;
}
}
