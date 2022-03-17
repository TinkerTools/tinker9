#include "integrator.h"
#include "md.h"
#include "tool/io.h"

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
   mdKinetic(temp);
   mdBussiThermostat(dt, temp);
}
}
