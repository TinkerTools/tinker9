#pragma once
#include "itgiVerlet.h"

namespace tinker {
class RespaIntegrator : public VerletIntegrator
{
public:
   ~RespaIntegrator();
   RespaIntegrator(ThermostatEnum, BarostatEnum);
   void kickoff() override;
};
}
