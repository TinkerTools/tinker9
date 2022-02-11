#pragma once
#include "itgtBasic.h"
#include "nhc.h"
#include <string>

namespace tinker {
class NhcThermostat : public BasicThermostat
{
protected:
   static constexpr int nhc_nsy = 3;
   int nnose, nhc_nc;
   double g0;
   double vnh[maxnose], qnh[maxnose];
   double* (*f_kin)();
   void (*scale_vel)(double);
   std::string name;

   void controlImpl(double timeStep);

public:
   NhcThermostat(int nhclen, int nc, double dfree, double* (*kin)(),
                 void (*scale)(double), std::string str);
   void printDetail(FILE*);
   void control1(time_prec time_prec) override;
   void control2(time_prec time_prec, bool) override;

   static double* atomicKinetic();
   static void scaleAtomicVelocity(double scale);
};
}
