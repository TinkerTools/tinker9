#include "box.h"
#include "energy.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include "mdpt.h"
#include "nose.h"
#include "random.h"
#include <cmath>
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/stodyn.hh>
#include <tinker/detail/units.hh>


namespace tinker {
namespace {
double sinh_id(double x)
{
   double y = std::fabs(x);
   if (y <= 1.0e-8)
      return 1.0;
   else
      return std::sinh(y) / y;
}


void lp_v5(time_prec dt, double R)
{
   const double D = 3.0;
   const double Nf = mdstuf::nfree;
   const double g = stodyn::friction;
   const double kbt = units::gasconst * bath::kelvin;
   const double odnf = 1.0 + D / Nf;


   const int nc = 2;
   const double h = dt / (2 * nc);
   const double h_2 = 0.5 * h;
   const double opgh2 = 1.0 + g * h_2;
   const double omgh2 = 1.0 - g * h_2;
   const double sdbar = std::sqrt(2.0 * kbt * g * dt / qbar) / (4 * nc);


   const double& qnh0 = qnh[0];
   double& gnh0 = gnh[0];
   double& vnh0 = vnh[0];


   const virial_prec tr_vir = vir[0] + vir[4] + vir[8];
   const double vol0 = volbox();
   double temp0;
   kinetic(temp0);


   double DelP;
   double eksum0 = eksum, eksum1;
   double velsc0 = 1.0, velsc1;
   for (int k = 0; k < nc; ++k) {
      eksum1 = eksum0;
      velsc1 = velsc0;


      // vnh 1/2
      gnh0 = (2 * eksum1 + qbar * vbar * vbar - Nf * kbt) / qnh0;
      vnh0 = vnh0 + gnh0 * h_2;


      // vbar 1/2
      DelP = (odnf * 2 * eksum1 - tr_vir);
      DelP = DelP - D * vol0 * bath::atmsph / units::prescon;
      gbar = DelP / qbar;
      vbar = omgh2 * vbar + gbar * h_2 + sdbar * R;


      // velocity
      double scal = std::exp(-h * (odnf * vbar + vnh0));
      velsc1 *= scal;
      eksum1 *= (scal * scal);


      // vbar 2/2
      DelP = (odnf * 2 * eksum1 - tr_vir);
      DelP = DelP - D * vol0 * bath::atmsph / units::prescon;
      gbar = DelP / qbar;
      vbar = (vbar + gbar * h_2 + sdbar * R) / opgh2;


      // vnh 2/2
      gnh0 = (2 * eksum1 + qbar * vbar * vbar - Nf * kbt) / qnh0;
      vnh0 = vnh0 + gnh0 * h_2;


      eksum0 = eksum1;
      velsc0 = velsc1;
   }


   darray::scale(g::q0, n, velsc0, vx);
   darray::scale(g::q0, n, velsc0, vy);
   darray::scale(g::q0, n, velsc0, vz);
   const double velsc2 = velsc0 * velsc0;
   eksum *= velsc2;
   for (int ii = 0; ii < 3; ++ii)
      for (int jj = 0; jj < 3; ++jj)
         ekin[ii][jj] *= velsc2;
}


void vv_lpiston_npt_v5(int istep, time_prec dt)
{
   int vers1 = rc_flag & calc::vmask;
   bool save = !(istep % inform::iwrite);
   if (!save)
      vers1 &= ~calc::energy;


   const time_prec dt_2 = 0.5 * dt;
   const double R = normal<double>();


   lp_v5(dt, R);
   propagate_velocity(dt_2, gx, gy, gz);


   // volume
   const double term = vbar * dt_2;
   const double expterm = std::exp(term);
   const double eterm2 = expterm * expterm;
   lvec1 *= eterm2;
   lvec2 *= eterm2;
   lvec3 *= eterm2;
   set_default_recip_box();


   // sinh(x)/x
   double poly = sinh_id(term);
   poly *= expterm * dt;
   propagate_pos_axbv(eterm2, poly);
   copy_pos_to_xyz(true);


   energy(vers1);


   propagate_velocity(dt_2, gx, gy, gz);
   lp_v5(dt, R);
}
}


void vv_lpiston_npt(int istep, time_prec dt)
{
   vv_lpiston_npt_v5(istep, dt);
}
}
