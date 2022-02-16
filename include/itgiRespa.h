#pragma once
#include "itgiBasic.h"

namespace tinker {
class RespaIntegrator : public BasicIntegrator
{
protected:
   void kickoff() override;

public:
   RespaIntegrator(ThermostatEnum te, BarostatEnum be);
   static void KickOff();
};
}
