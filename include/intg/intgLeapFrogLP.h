#pragma once
#include "intgBasic.h"

namespace tinker {
class LeapFrogLPIntegrator : public BasicIntegrator
{
public:
   ~LeapFrogLPIntegrator();
   LeapFrogLPIntegrator();
   void printDetail(FILE*) override;
   void kickoff() override;
   void dynamic(int, time_prec) override;
};
}
