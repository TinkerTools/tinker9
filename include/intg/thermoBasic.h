#pragma once
#include "mdprec.h"
#include <cstdio>

namespace tinker {
enum class ThermostatEnum;

class BasicThermostat
{
protected:
   ThermostatEnum m_thermoEnum;
   void printBasic(FILE*);
   BasicThermostat(ThermostatEnum);

public:
   BasicThermostat();
   virtual ~BasicThermostat();
   virtual void printDetail(FILE*);
   virtual void control1(time_prec timeStep) {}
   virtual void control2(time_prec timeStep, bool save);
};
}
