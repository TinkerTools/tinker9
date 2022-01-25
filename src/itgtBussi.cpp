#include "itgtBussi.h"
#include "itgEnum.h"
#include "mdpt.h"

namespace tinker {
BussiThermostat::BussiThermostat()
   : BasicThermostat(ThermostatEnum::Bussi)
{}

void BussiThermostat::control2(time_prec dt, bool)
{
   double temp;
   kinetic(temp);
   bussi_thermostat(dt, temp);
}
}
