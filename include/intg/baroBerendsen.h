#pragma once
#include "baroBasic.h"

namespace tinker {
class BerendsenBarostat : public BasicBarostat
{
public:
   void control2(time_prec timeStep) override;
};
}
