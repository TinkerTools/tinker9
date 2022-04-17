#include "md/integrator.h"
#include "md/misc.h"
#include "tool/ioprint.h"
#include <tinker/detail/mdstuf.hh>

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

BasicThermostat* BasicThermostat::create(ThermostatEnum te)
{
   BasicThermostat* t = nullptr;
   switch (te) {
   case ThermostatEnum::BUSSI:
      t = new BussiThermostat;
      break;
   case ThermostatEnum::NHC:
      t = new NhcDevice(5, 5, static_cast<double>(mdstuf::nfree), //
         NhcDevice::kineticAtomic,                                //
         NhcDevice::scaleVelocityAtomic,                          //
         std::string("NHC"));
      break;
   default:
      t = new BasicThermostat;
      break;
   }
   return t;
}
}
