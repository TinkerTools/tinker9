#pragma once
#include "itgtBasic.h"

namespace tinker {
class BussiThermostat : public BasicThermostat
{
public:
   void control2(double timeStep, bool) override;
};
}
