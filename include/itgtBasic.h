#pragma once
#include "itgEnum.h"
#include "mdprec.h"
#include <cstdio>

namespace tinker {
class BasicThermostat : public IntegratorStaticData
{
protected:
   void printBasic(FILE*);

public:
   BasicThermostat();
   virtual ~BasicThermostat();
   virtual void printDetail(FILE*);
   virtual void control1(time_prec timeStep);
   virtual void control2(time_prec timeStep, bool save);
};
}
