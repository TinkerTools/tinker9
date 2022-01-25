#pragma once
#include "itgbBasic.h"

namespace tinker {
class MonteCarloBarostat : public BasicBarostat
{
public:
   ~MonteCarloBarostat();
   MonteCarloBarostat();
   void control4(time_prec) override;
   bool ifApply(int istep) override;
};
}
