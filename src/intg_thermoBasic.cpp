#include "intg/thermoBasic.h"
#include "macro.h"
#include "tool/io_print.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>

namespace tinker {
void BasicThermostat::printBasic(FILE* o)
{
   if (not inform::verbose)
      return;

   print(o, " Temperature      : %12.1lf\n", bath::kelvin);
   print(o, " Tau-Temperature  : %12.1lf\n", bath::tautemp);
   print(o, "\n");
}

BasicThermostat::BasicThermostat() {}

BasicThermostat::~BasicThermostat() {}

void BasicThermostat::printDetail(FILE* o)
{
   printBasic(o);
}
}
