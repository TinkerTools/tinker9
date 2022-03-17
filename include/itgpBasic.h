#pragma once
#include "itgEnum.h"
#include "md.h"

namespace tinker {
class BasicPropagator : public IntegratorStaticData
{
public:
   bool ifSave(int istep) const;
   BasicPropagator();
   virtual ~BasicPropagator();

   virtual void updateVelocity1(time_prec t);
   virtual void updateVelocity2(time_prec t);
   virtual void updateVelocityR0(time_prec t);
   virtual void updateVelocityR1(time_prec t, int nrespa);
   virtual void updateVelocityR2(time_prec t, int nrespa);
   virtual void updatePosition(time_prec t);

   virtual void rattleSave();
   virtual void rattle(time_prec dt);
   virtual void rattle2(time_prec dt, bool useVirial);
   static bool useRattle();
};

typedef BasicPropagator VerletPropagator;
}
