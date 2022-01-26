#include "itgEnum.h"
#include "itgbBasic.h"
#include "itgbBerendsen.h"
#include "itgbMonteCarlo.h"
#include "itgtBasic.h"
#include "itgtBussi.h"
#include "itgtNhc96.h"
#include "mdpq.h"
#include <tinker/detail/mdstuf.hh>

namespace tinker {
BasicThermostat* create(ThermostatEnum te)
{
   // TODO check nfree
   if (te == ThermostatEnum::Bussi)
      return new BussiThermostat;
   else if (te == ThermostatEnum::Nhc96)
      return new Nhc96Thermostat(5, 5, static_cast<double>(mdstuf::nfree),
                                 Nhc96Thermostat::atomicKinetic,
                                 Nhc96Thermostat::scaleAtomicVelocity,
                                 std::string("NHC"));
   else
      return new BasicThermostat;
}

BasicBarostat* create(BarostatEnum be)
{
   if (be == BarostatEnum::Berendsen)
      return new BerendsenBarostat;
   else if (be == BarostatEnum::MonteCarlo)
      return new MonteCarloBarostat;
   else
      return new BasicBarostat;
}
}
