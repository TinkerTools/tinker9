#include "md/integrator.h"
#include "md/misc.h"
#include "tool/ioprint.h"

namespace tinker {
void BussiThermostat::printDetail(FILE* o)
{
   print(o, "\n");
   print(o, " %s\n", "Bussi Thermostat");
   printBasic(o);
}

void BussiThermostat::control2(time_prec dt, bool)
{
   double temp;
   kinetic(temp);
   bussiThermostat(dt, temp);
}
}
