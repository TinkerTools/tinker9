#pragma once
#include "itgbBasic.h"

namespace tinker {
class BerendsenBarostat : public BasicBarostat
{
public:
   BarostatEnum getBarostatEnum() const override;
   void control2(time_prec timeStep) override;
};
}
