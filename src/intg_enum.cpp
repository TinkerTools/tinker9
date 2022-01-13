#include "intg/baroBasic.h"
#include "intg/baroBerendsen.h"
#include "intg/baroMonteCarlo.h"
#include "intg/enum.h"
#include "intg/thermoBasic.h"

namespace tinker {
BasicThermostat* create(ThermostatEnum te)
{
   if (te == ThermostatEnum::Null)
      return new BasicThermostat;
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
