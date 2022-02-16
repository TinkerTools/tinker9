#pragma once
#include "itgiBasic.h"

namespace tinker {
class VerletIntegrator : public BasicIntegrator
{
protected:
   void kickoff() override;

public:
   VerletIntegrator(ThermostatEnum te, BarostatEnum be);
   static void KickOff();
};
}
