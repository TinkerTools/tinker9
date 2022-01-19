#pragma once
#include "thermoBasic.h"

namespace tinker {
class BussiThermostat : public BasicThermostat
{
public:
   BussiThermostat();
   void control2(double timeStep, bool) override;
};
}
