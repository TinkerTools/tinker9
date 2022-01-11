#pragma once
#include <cstdio>

namespace tinker {
class BasicThermostat
{
protected:
   void printBasic(FILE*);

public:
   BasicThermostat();
   virtual ~BasicThermostat();
   virtual void printDetail(FILE*);
   virtual void control(double timeStep);
   virtual void control1(double timeStep);
   virtual void control2(double timeStep);
};
}
