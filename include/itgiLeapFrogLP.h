#pragma once
#include "itgiBasic.h"

namespace tinker {
class LeapFrogLPIntegrator : public BasicIntegrator
{
protected:
   const char* name() const override;
   void kickoff() override;

public:
   LeapFrogLPIntegrator();
   ~LeapFrogLPIntegrator();
   void dynamic(int, time_prec) override;
};
}
