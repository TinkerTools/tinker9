#pragma once
#include "itgiBasic.h"

namespace tinker {
class VerletIntegrator : public BasicIntegrator
{
protected:
   const char* name() const override;
   void kickoff() override;

public:
   VerletIntegrator(ThermostatEnum te, BarostatEnum be);
   static void KickOff();
};
}
