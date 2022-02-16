#include "itgiNhc06.h"
#include "itgiRespa.h"
#include "itgiVerlet.h"
#include "itgpLogV.h"
#include "tool/error.h"
#include <tinker/detail/bath.hh>

namespace tinker {
void Nhc06Integrator::kickoff()
{
   if (m_isNRespa1)
      VerletIntegrator::KickOff();
   else
      RespaIntegrator::KickOff();
}

Nhc06Integrator::Nhc06Integrator(bool isNRespa1)
   : BasicIntegrator(PropagatorEnum::Verlet, ThermostatEnum::m_Nhc2006,
                     BarostatEnum::Nhc2006)
   , m_isNRespa1(isNRespa1)
{
   if (BasicPropagator::useRattle())
      TINKER_THROW("Constraints under NH-NPT require the ROLL algorithm.");

   if (bath::anisotrop)
      TINKER_THROW("Cannot use ANISO-PRESSURE in Nhc06Integrator.");

   bool isAtomic = true;
   bool isAniso = false;
   bool isPedantic = true;

   delete m_prop;
   m_prop = new LogVPropagator(isNRespa1, isAtomic, isAniso, isPedantic);

   this->kickoff();
}
}
