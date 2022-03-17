#pragma once
#include "itgtBasic.h"
#include "nhc.h"
#include <string>

namespace tinker {
/**
 * \ingroup mdpt
 * \brief Applies a velocity correction as needed for the Nose-Hoover Chains
 * at the half time step.
 *
 * Literature reference:
 *    - <a href="https://doi.org/10.1080/00268979600100761">
 *    G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein,
 *    "Explicit Reversible Integrators for Extended Systems Dynamics",
 *    Molecular Physics, 87, 1117-1157 (1996).
 *    </a>
 */
class NhcDevice : public BasicThermostat
{
protected:
   static constexpr int nhc_nsy = 3;
   int nnose, nhc_nc;
   double g0;
   double vnh[maxnose], qnh[maxnose];
   double (*f_kin)();
   void (*scale_vel)(double);
   std::string name;

   void controlImpl(double timeStep);

public:
   NhcDevice(int nhclen, int nc, double dfree, //
      double (*kin)(), void (*scale)(double), std::string str);
   void printDetail(FILE*);
   void control1(time_prec time_prec) override;
   void control2(time_prec time_prec, bool) override;

   static double kineticAtomic();
   static void scaleVelocityAtomic(double scale);
};
}
