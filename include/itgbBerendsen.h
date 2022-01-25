#pragma once
#include "itgbBasic.h"

namespace tinker {
class BerendsenBarostat : public BasicBarostat
{
public:
   BerendsenBarostat();
   void control2(time_prec timeStep) override;
};
}
