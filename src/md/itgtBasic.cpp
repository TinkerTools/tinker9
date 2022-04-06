#include "md/integrator.h"
#include "md/pt.h"
#include "tool/ioprint.h"
#include <tinker/detail/bath.hh>

namespace tinker {
void BasicThermostat::printBasic(FILE* o)
{
   print(o, " Temperature        %12.1lf Kelvin\n", bath::kelvin);
   print(o, " Tau-Temperature    %12.1lf ps\n", bath::tautemp);
}

BasicThermostat::BasicThermostat() {}

BasicThermostat::~BasicThermostat() {}

void BasicThermostat::printDetail(FILE* o) {}

void BasicThermostat::control1(time_prec) {}

void BasicThermostat::control2(time_prec, bool save)
{
   if (save) {
      T_prec temp;
      kinetic(temp);
   }
}
}
