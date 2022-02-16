#pragma once
#include "itgiBasic.h"

namespace tinker {
class LeapFrogLPIntegrator : public BasicIntegrator
{
protected:
   void kickoff() override;

public:
   LeapFrogLPIntegrator();
   ~LeapFrogLPIntegrator();
   void printDetail(FILE*) override;
   void dynamic(int, time_prec) override;
};
}
