#pragma once
#include "itgiBasic.h"

namespace tinker {
class Nhc96Integrator : public BasicIntegrator
{
protected:
   void kickoff() override;

public:
   Nhc96Integrator();
   void printDetail(FILE*) override;
   void dynamic(int, time_prec) override;
};
}
