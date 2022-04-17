#include "ff/energy.h"
#include "md/integrator.h"
#include "md/pq.h"
#include "md/rattle.h"
#include "tool/argkey.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>

namespace tinker {
bool IntegratorStaticData::applyBaro = false;
bool IntegratorStaticData::atomic = true;
bool IntegratorStaticData::aniso = false;
bool IntegratorStaticData::semiiso = false;
int IntegratorStaticData::nrespa = 0;
double IntegratorStaticData::dofP = 1.0;
int IntegratorStaticData::anisoArrayLength = 0;
const int IntegratorStaticData::anisoArray[6][2] = {{2, 2}, {1, 1}, {0, 0}, {0, 2}, {0, 1}, {1, 2}};
}

namespace tinker {
BasicPropagator::BasicPropagator()
{
   nrespa = 1;

   atomic = not useRattle();
   aniso = bath::anisotrop;
   getKV("SEMIISO-PRESSURE", semiiso, false);
   if (semiiso)
      aniso = true;
}

BasicPropagator::~BasicPropagator() {}

bool BasicPropagator::ifSave(int istep) const
{
   return (istep % inform::iwrite) == 0;
}

void BasicPropagator::pos(time_prec t)
{
   mdPos(t);
}

void BasicPropagator::vel1(time_prec t)
{
   mdVel(t, gx, gy, gz);
}

void BasicPropagator::vel2(time_prec t)
{
   this->vel1(t);
}

void BasicPropagator::rattleSave()
{
   if (not useRattle())
      return;

   darray::copy(g::q0, n, rattle_xold, xpos);
   darray::copy(g::q0, n, rattle_yold, ypos);
   darray::copy(g::q0, n, rattle_zold, zpos);
}

void BasicPropagator::rattle(time_prec timeStep)
{
   if (not useRattle())
      return;

   tinker::rattle(timeStep, rattle_xold, rattle_yold, rattle_zold);
}

void BasicPropagator::rattle2(time_prec timeStep, bool useVirial)
{
   if (not useRattle())
      return;

   tinker::rattle2(timeStep, useVirial);
}

BasicPropagator* BasicPropagator::create(PropagatorEnum pe)
{
   BasicPropagator* p = nullptr;
   switch (pe) {
   case PropagatorEnum::RESPA:
      p = new RespaDevice;
      break;
   default:
      p = new BasicPropagator;
      break;
   }
   return p;
}
}
