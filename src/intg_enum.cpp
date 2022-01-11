#include "intg/baroBasic.h"
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
   if (be == BarostatEnum::Null)
      return new BasicBarostat;
   else
      return new BasicBarostat;
}
}
