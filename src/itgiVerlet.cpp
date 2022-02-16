#include "itgiVerlet.h"
#include "energy.h"
#include "mdcalc.h"
#include "mdpq.h"

namespace tinker {
void VerletIntegrator::kickoff()
{
   VerletIntegrator::KickOff();
}

VerletIntegrator::VerletIntegrator(ThermostatEnum te, BarostatEnum be)
   : BasicIntegrator(PropagatorEnum::Verlet, te, be)
{
   this->kickoff();
}

void VerletIntegrator::KickOff()
{
   energy((calc::grad | calc::virial) & rc_flag);
}
}
