#include "intg/enum.h"
#include "intg/thermoBasic.h"
#include "mdpt.h"
#include "tool/io_print.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>

namespace tinker {
void BasicThermostat::printBasic(FILE* o)
{
   if (not inform::verbose)
      return;

   print(o, " Temperature      : %12.1lf Kelvin\n", bath::kelvin);
   print(o, " Tau-Temperature  : %12.1lf ps\n", bath::tautemp);
   print(o, "\n");
}

BasicThermostat::BasicThermostat(ThermostatEnum te)
   : m_thermoEnum(te)
{}

BasicThermostat::BasicThermostat()
   : m_thermoEnum(ThermostatEnum::Null)
{}

BasicThermostat::~BasicThermostat() {}

void BasicThermostat::printDetail(FILE* o)
{
   printBasic(o);
}

void BasicThermostat::control2(time_prec, bool save)
{
   if (save) {
      T_prec temp;
      kinetic(temp);
   }
}
}
