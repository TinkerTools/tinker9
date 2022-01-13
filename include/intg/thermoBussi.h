#pragma once
#include "thermoBasic.h"

namespace tinker {
class BussiThermostat : public BasicThermostat
{
public:
   void control2(double timeStep) override;
};
}
