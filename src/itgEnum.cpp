#include "integrator.h"

namespace tinker {
BasicPropagator* create(PropagatorEnum pe)
{
   BasicPropagator* p = nullptr;
   switch (pe) {
   case PropagatorEnum::Respa:
      p = new RespaDevice;
      break;
   default:
      p = new BasicPropagator;
      break;
   }
   return p;
}
}

#include "md.h"
#include <tinker/detail/mdstuf.hh>
namespace tinker {
BasicThermostat* create(ThermostatEnum te)
{
   BasicThermostat* t = nullptr;
   switch (te) {
   case ThermostatEnum::Bussi:
      t = new BussiThermostat;
      break;
   case ThermostatEnum::Nhc:
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
   case BarostatEnum::Berendsen:
      b = new BerendsenBarostat;
      break;
   case BarostatEnum::LP2022:
      b = new LP22Barostat;
      break;
   case BarostatEnum::MonteCarlo:
      b = new MonteCarloBarostat;
      break;
   case BarostatEnum::Nhc2006:
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
const int IntegratorStaticData::anisoArray[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 2}, {0, 1}, {1, 2}};
}
