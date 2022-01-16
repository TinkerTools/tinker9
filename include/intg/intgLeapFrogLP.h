#pragma once
#include "intgBasic.h"

namespace tinker {
class LeapFrogLPIntegrator : public BasicIntegrator
{
public:
   ~LeapFrogLPIntegrator();
   LeapFrogLPIntegrator();
   void kickoff() override;
   void dynamic(int, time_prec) override;
};
}
