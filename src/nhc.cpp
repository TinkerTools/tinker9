#include "nhc.h"
#include <cassert>
#include <cmath>
#include <tinker/detail/units.hh>

namespace tinker {
void nhc_isot_96(time_prec dt, int nnose, double* vnh, const double* qnh, double g0,
   double* (*f_kin)(), void (*scale_vel)(double))
{
   constexpr int nc = 5;
   constexpr int ns = 3;
   // w[0] = 1/(2 - 2**(1/3))
   // w[1] = 1 - w[0] - w[2]
   // w[2] = w[0]
   constexpr double w[3] = {1.351207191959657634047687808971460826921999376217144828328,
      -1.70241438391931526809537561794292165384399875243428965665,
      1.351207191959657634047687808971460826921999376217144828328};

   assert(nnose <= maxnose);
   double gnh[maxnose];

   const double kbt = units::gasconst * bath::kelvin;
   double& ref_eksum = *f_kin();

   const time_prec dtc = dt / nc;
   double eksum0 = ref_eksum;
   double velsc0 = 1.0;
   for (int k = 0; k < nc; ++k) {
      for (int j = 0; j < ns; ++j) {
         const time_prec dts = w[j] * dtc;
         const time_prec dt2 = 0.5 * dts;
         const time_prec dt4 = 0.25 * dts;
         const time_prec dt8 = 0.125 * dts;

         for (int i = nnose - 1; i > -1; --i) {
            if (i == 0)
               gnh[i] = (2 * eksum0 - g0 * kbt) / qnh[i];
            else
               gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - kbt) / qnh[i];

            if (i == nnose - 1)
               vnh[i] += gnh[i] * dt4;
            else {
               double exptm = std::exp(-vnh[i + 1] * dt8);
               vnh[i] = (vnh[i] * exptm + gnh[i] * dt4) * exptm;
            }
         }

         double scal;
         scal = std::exp(-dt2 * vnh[0]);
         velsc0 *= scal;
         eksum0 *= (scal * scal);

         for (int i = 0; i < nnose; ++i) {
            if (i == 0)
               gnh[i] = (2 * eksum0 - g0 * kbt) / qnh[i];
            else
               gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - kbt) / qnh[i];

            if (i == nnose - 1)
               vnh[i] += gnh[i] * dt4;
            else {
               double exptm = std::exp(-vnh[i + 1] * dt8);
               vnh[i] = (vnh[i] * exptm + gnh[i] * dt4) * exptm;
            }
         }
      }
   }

   ref_eksum = eksum0;
   scale_vel(velsc0);
}
}
