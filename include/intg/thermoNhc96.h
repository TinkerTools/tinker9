#pragma once
#include "nhc.h"
#include "thermoBasic.h"
#include <string>

namespace tinker {
class Nhc96Thermostat : public BasicThermostat
{
protected:
   static constexpr int nhc_nsy = 3;
   int nnose, nhc_nc;
   double g0;
   double vnh[maxnose], qnh[maxnose];
   double* (*f_kin)();
   void (*scale_vel)(double);
   void controlImpl(double timeStep);
   std::string name;

public:
   Nhc96Thermostat(int nhclen, int nc, double dfree, double* (*kin)(),
                   void (*scale)(double), std::string str);
   void printDetail(FILE*);
   void control1(time_prec time_prec) override;
   void control2(time_prec time_prec, bool) override;
};
}
