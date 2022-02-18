#pragma once
#include "itgiBasic.h"

namespace tinker {
class Nhc96Integrator : public BasicIntegrator
{
protected:
   const char* name() const override;
   void kickoff() override;

public:
   Nhc96Integrator();
   void dynamic(int, time_prec) override;
};
}
