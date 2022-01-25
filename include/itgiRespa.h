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
};
}
