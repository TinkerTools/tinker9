#pragma once
#include "itgpRespa.h"

namespace tinker {
class LogVPropagator : public BasicPropagator
{
protected:
   typedef BasicPropagator base_t;

   RespaPropagator* m_respa;
   bool m_atomic;   // atomic vs. molecular
   bool m_aniso;    // anisotropic cell fluctuation
   bool m_pedantic; // pedantic velocity update

public:
   LogVPropagator(bool isNRespa1, bool isAtomic, bool isAniso, bool isPedantic);
   ~LogVPropagator();

   void updatePosition(time_prec t) override;

   void updateVelocity0(time_prec t) override;
   void updateVelocity1(time_prec t) override;
   void updateVelocity2(time_prec t) override;

   void updateVelocityR1(time_prec tfast, time_prec t) override;
   void updateVelocityR2(time_prec tfast, time_prec t) override;
};
}
