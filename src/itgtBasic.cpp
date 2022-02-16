#include "itgtBasic.h"
#include "mdpt.h"
#include "tool/io_print.h"
#include <tinker/detail/bath.hh>

namespace tinker {
void BasicThermostat::printBasic(FILE* o)
{
   print(o, "\n");
   print(o, " Temperature      : %12.1lf Kelvin\n", bath::kelvin);
   print(o, " Tau-Temperature  : %12.1lf ps\n", bath::tautemp);
}

BasicThermostat::BasicThermostat() {}

BasicThermostat::~BasicThermostat() {}

void BasicThermostat::printDetail(FILE* o)
{
   printBasic(o);
}

void BasicThermostat::control1(time_prec) {}

void BasicThermostat::control2(time_prec, bool save)
{
   if (save) {
      T_prec temp;
      kinetic(temp);
   }
}
}
