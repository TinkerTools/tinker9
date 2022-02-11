#pragma once
#include "itgiVerlet.h"

namespace tinker {
class RespaIntegrator : public VerletIntegrator
{
protected:
   void kickoff() override;

public:
   ~RespaIntegrator();
   RespaIntegrator(ThermostatEnum, BarostatEnum);
   RespaIntegrator(PropagatorEnum pe, ThermostatEnum te, BarostatEnum be) : VerletIntegrator(pe,te,be) {}
};
}
