#include "itgtBussi.h"
#include "mdpt.h"
#include "tool/io_print.h"

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
   bussi_thermostat(dt, temp);
}
}
