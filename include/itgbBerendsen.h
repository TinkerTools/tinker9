#pragma once
#include "itgbBasic.h"

namespace tinker {
class BerendsenBarostat : public BasicBarostat
{
public:
   void printDetail(FILE*) override;
   BarostatEnum getBarostatEnum() const override;
   void control2(time_prec timeStep) override;
};
}
