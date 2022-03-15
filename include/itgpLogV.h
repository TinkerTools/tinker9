#pragma once
#include "itgpRespa.h"

namespace tinker {
class LogVPropagator : public BasicPropagator
{
protected:
   typedef BasicPropagator base_t;

   RespaPropagator* m_respa;

   void updateVelocityImpl(time_prec t, int idx, int nrespa);

public:
   LogVPropagator(bool isNRespa1);
   ~LogVPropagator();

   void updatePosition(time_prec t) override;
   void updateVelocity1(time_prec t) override;
   void updateVelocity2(time_prec t) override;

   void updateVelocityR0(time_prec t) override;
   void updateVelocityR1(time_prec t, int nrespa) override;
   void updateVelocityR2(time_prec t, int nrespa) override;
};
}
