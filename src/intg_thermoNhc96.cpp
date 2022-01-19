#include "intg/enum.h"
#include "intg/thermoNhc96.h"
#include <algorithm>
#include <cmath>
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/units.hh>

namespace tinker {
void Nhc96Thermostat::controlImpl(time_prec dt)
{
   const int nc = nhc_nc;
   constexpr int ns = nhc_nsy;
   static_assert(ns == 3, "");
   constexpr double w[3] = {
      1.351207191959657634047687808971460826921999376217144828328,
      -1.70241438391931526809537561794292165384399875243428965665,
      1.351207191959657634047687808971460826921999376217144828328};

   double gnh[maxnose];
   double kbt = units::gasconst * bath::kelvin;
   double& ref_eksum = *f_kin();

   double dtc = dt / nc;
   double eksum0 = ref_eksum;
   double velsc0 = 1.0;
   for (int k = 0; k < nc; ++k) {
      for (int j = 0; j < ns; ++j) {
         double dts = w[j] * dtc;
         double dt2 = 0.5 * dts;
         double dt4 = 0.25 * dts;
         double dt8 = 0.125 * dts;

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

Nhc96Thermostat::Nhc96Thermostat(int nhclen, int nc, double dfree,
                                 double* (*kin)(), void (*scale)(double),
                                 std::string str)
   : BasicThermostat(ThermostatEnum::Nhc96)
   , nnose(nhclen)
   , nhc_nc(nc)
   , g0(dfree)
   , f_kin(kin)
   , scale_vel(scale)
   , name(str)
{
   nnose = std::max(1, nnose);
   nnose = std::min(nnose, maxnose);

   nhc_nc = std::max(1, nhc_nc);
   nhc_nc = std::min(nhc_nc, 5);

   // default vnh and qnh
   double kt = units::gasconst * bath::kelvin;
   for (int i = 0; i < maxnose; ++i) {
      qnh[i] = kt * bath::tautemp * bath::tautemp;
      vnh[i] = 0;
   }
   qnh[0] *= g0;
}

void Nhc96Thermostat::printDetail(FILE* o)
{
   if (not inform::verbose)
      return;

   fprintf(o, " %s\n", name.c_str());
   fprintf(o, " NHC N            : %12d\n", nnose);
   fprintf(o, " NHC NC           : %12d\n", nhc_nc);
   fprintf(o, " NHC NSY          : %12d\n", nhc_nsy);
   for (int i = 0; i < nnose; ++i) {
      fprintf(o, " NHC %2d Mass      : %12.4lf\n", i + 1, qnh[i]);
      fprintf(o, " NHC %2d Vel       : %12.4lf\n", i + 1, vnh[i]);
   }
   printBasic(o);
}

void Nhc96Thermostat::control1(time_prec dt)
{
   this->controlImpl(dt);
}

void Nhc96Thermostat::control2(time_prec dt, bool)
{
   this->controlImpl(dt);
}
}
