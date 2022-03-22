#include "md/integrator.h"
#include "tool/error.h"
#include <tinker/detail/bath.hh>

namespace tinker {
const char* Nhc06Integrator::name() const
{
   return "NHC2006";
}

void Nhc06Integrator::kickoff()
{
   if (m_isNRespa1)
      VerletIntegrator::KickOff();
   else
      RespaIntegrator::KickOff();
}

Nhc06Integrator::Nhc06Integrator(bool isNRespa1)
   : BasicIntegrator(PropagatorEnum::Verlet, ThermostatEnum::m_Nhc2006, BarostatEnum::Nhc2006)
   , m_isNRespa1(isNRespa1)
{
   if (useRattle())
      TINKER_THROW("Constraints under NH-NPT require the ROLL algorithm.");

   if (bath::anisotrop)
      TINKER_THROW("Cannot use ANISO-PRESSURE in Nhc06Integrator.");

   delete m_prop;
   m_prop = new LogVDevice(isNRespa1);

   this->kickoff();
}
}
