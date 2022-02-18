#pragma once
#include "itgiBasic.h"

namespace tinker {
class RespaIntegrator : public BasicIntegrator
{
protected:
   const char* name() const override;
   void kickoff() override;

public:
   RespaIntegrator(ThermostatEnum te, BarostatEnum be);
   static void KickOff();
};
}
