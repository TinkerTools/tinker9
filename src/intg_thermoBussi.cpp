#include "intg/thermoBussi.h"
#include "mdpt.h"

namespace tinker {
void BussiThermostat::control2(time_prec dt)
{
   double temp;
   kinetic(temp);
   bussi_thermostat(dt, temp);
}
}
