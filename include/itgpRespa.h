#pragma once
#include "itgpBasic.h"

namespace tinker {
class RespaPropagator : public BasicPropagator
{
public:
   RespaPropagator();
   ~RespaPropagator();

   void updateVelocityR1(time_prec t, int nrespa) override;
   void updateVelocityR2(time_prec t, int nrespa) override;
};
}
