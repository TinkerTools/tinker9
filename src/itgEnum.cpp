#include "itgEnum.h"

#include "itgpBasic.h"
#include "itgpRespa.h"
namespace tinker {
BasicPropagator* create(PropagatorEnum pe)
{
   BasicPropagator* p = nullptr;
   switch (pe) {
   case PropagatorEnum::Respa:
      p = new RespaPropagator;
      break;
   default:
      p = new BasicPropagator;
      break;
   }
   return p;
}
}

#include "itgtBasic.h"
#include "itgtBussi.h"
#include "itgtNhc.h"
#include "mdpq.h"
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
      t = new NhcThermostat(5, 5, static_cast<double>(mdstuf::nfree), //
                            NhcThermostat::kineticAtomic,             //
                            NhcThermostat::scaleVelocityAtomic,       //
                            std::string("NHC"));
      break;
   default:
      t = new BasicThermostat;
      break;
   }
   return t;
}
}

#include "itgbBasic.h"
#include "itgbBerendsen.h"
#include "itgbLP22.h"
#include "itgbMonteCarlo.h"
#include "itgbNhc06.h"
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
double IntegratorStaticData::dofP = 1.0;
int IntegratorStaticData::anisoArrayLength = 0;
const int IntegratorStaticData::anisoArray[6][2] = {{0, 0}, {1, 1}, {2, 2},
                                                    {0, 2}, {0, 1}, {1, 2}};

void __PlaceHolderMessage(const char* s)
{
   printf("\x1b[31m%s\033[0m\n", s);
}
}
