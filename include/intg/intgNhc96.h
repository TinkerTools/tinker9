#pragma once
#include "intgBasic.h"

namespace tinker {
class Nhc96Integrator : public BasicIntegrator
{
public:
   void printDetail(FILE*) override;
   void kickoff() override;
   void dynamic(int, time_prec) override;
};
}
