#pragma once
#include "itgbBasic.h"

namespace tinker {
class MonteCarloBarostat : public BasicBarostat
{
public:
   ~MonteCarloBarostat();
   MonteCarloBarostat();
   void printDetail(FILE*) override;
   BarostatEnum getBarostatEnum() const override;
   void control4(time_prec) override;
   bool ifApply(int istep) override;
};
}
