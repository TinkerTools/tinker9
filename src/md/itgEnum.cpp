#include "md/integrator.h"

namespace tinker {
BasicPropagator* create(PropagatorEnum pe)
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

#include <tinker/detail/mdstuf.hh>
namespace tinker {
BasicThermostat* create(ThermostatEnum te)
{
   BasicThermostat* t = nullptr;
   switch (te) {
   case ThermostatEnum::BUSSI:
      t = new BussiThermostat;
      break;
   case ThermostatEnum::NHC:
      t = new NhcDevice(5, 5, static_cast<double>(mdstuf::nfree), //
         NhcDevice::kineticAtomic,                                //
         NhcDevice::scaleVelocityAtomic,                          //
         std::string("NHC"));
      break;
   default:
      t = new BasicThermostat;
      break;
   }
   return t;
}
}

namespace tinker {
BasicBarostat* create(BarostatEnum be)
{
   BasicBarostat* b = nullptr;
   switch (be) {
   case BarostatEnum::BERENDSEN:
      b = new BerendsenBarostat;
      break;
   case BarostatEnum::LP2022:
      b = new LP22Barostat;
      break;
   case BarostatEnum::MONTECARLO:
      b = new MonteCarloBarostat;
      break;
   case BarostatEnum::NHC2006:
      b = new Nhc06Barostat;
      break;
   default:
      b = new BasicBarostat;
      break;
   }
   return b;
}
}

namespace tinker {
int IntegratorStaticData::nrespa = 0;
bool IntegratorStaticData::applyBaro = false;
bool IntegratorStaticData::atomic = true;
bool IntegratorStaticData::aniso = false;
bool IntegratorStaticData::semiiso = false;
double IntegratorStaticData::dofP = 1.0;
int IntegratorStaticData::anisoArrayLength = 0;
const int IntegratorStaticData::anisoArray[6][2] = {{2, 2}, {1, 1}, {0, 0}, {0, 2}, {0, 1}, {1, 2}};
}
