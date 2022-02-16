#include "itgtBussi.h"
#include "mdpt.h"

namespace tinker {
void BussiThermostat::control2(time_prec dt, bool)
{
   double temp;
   kinetic(temp);
   bussi_thermostat(dt, temp);
}
}
