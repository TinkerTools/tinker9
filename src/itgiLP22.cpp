#include "itgiLP22.h"
#include "itgiRespa.h"
#include "itgiVerlet.h"
#include "itgpLogV.h"
#include "tinker_rt.h"
#include <tinker/detail/bath.hh>

namespace tinker {
void LP22Integrator::kickoff()
{
   if (m_isNRespa1)
      VerletIntegrator::KickOff();
   else
      RespaIntegrator::KickOff();
}

LP22Integrator::LP22Integrator(bool isNRespa1)
   : BasicIntegrator(PropagatorEnum::Verlet, ThermostatEnum::m_LP2022,
                     BarostatEnum::LP2022)
   , m_isNRespa1(isNRespa1)
{
   bool aniso = bath::anisotrop;
   bool pedantic;
   get_kbool("PEDANTIC", pedantic, false);

   delete m_prop;
   m_prop = new LogVPropagator(isNRespa1, BasicPropagator::useRattle(), aniso,
                               pedantic);

   this->kickoff();
}
}
