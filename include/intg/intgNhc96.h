#pragma once
#include "intgBasic.h"

namespace tinker {
class Nhc96Integrator : public BasicIntegrator
{
public:
   void kickoff() override;
   void dynamic(int, time_prec) override;
};
}
