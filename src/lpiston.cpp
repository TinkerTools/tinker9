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
void lp(time_prec dt, double R)
{
   constexpr int D = 3;
   constexpr int nc = 5;
   constexpr int ns = 3;
   // w[0] = 1/(2 - 2**(1/3))
   // w[1] = 1 - w[0] - w[2]
   // w[2] = w[0]
   constexpr double w[3] = {
      1.351207191959657634047687808971460826921999376217144828328,
      -1.70241438391931526809537561794292165384399875243428965665,
      1.351207191959657634047687808971460826921999376217144828328};


   // trace of virial: xx+yy+zz
   const virial_prec tr_vir = vir[0] + vir[4] + vir[8];
   T_prec temp;
   kinetic(temp);
   double eksum0 = eksum;   // new kinetic energy
   double vbox0 = volbox(); // new box volume
   double vscal0 = 1.0;     // scalar for velocity
   double lensc0 = 1.0;     // scalar for box length


   for (int k = 0; k < nc; ++k) {
      for (int j = 0; j < ns; ++j) {
         double eksum1 = eksum0;
         double vbox1 = vbox0;
         double vscal1 = vscal0;
         double lensc1 = lensc0;


         const time_prec h = w[j] * dt / (2 * nc);
         const time_prec h_2 = h * 0.5;
         const double gh_2 = stodyn::friction * h_2;
         const double kbt = units::gasconst * bath::kelvin;
         const double sd = std::sqrt(kbt * gh_2 / qbar); // standard deviation


         const double odnf = 1.0 + D / mdstuf::nfree;
         const double opgh2 = 1.0 + gh_2;
         const double omgh2 = 1.0 - gh_2;


         // units::prescon ~ 6.86*10^4
         // 1 kcal/mol/Ang**3 = prescon atm
         double DelP, DelP_m;


         // BBK vbar 1/2
         DelP = odnf * 2 * eksum1 - tr_vir;
         DelP = DelP - D * vbox1 * bath::atmsph / units::prescon;
         DelP_m = DelP / qbar;
         vbar = DelP_m * h_2 + omgh2 * vbar + sd * R;


         // v and volume
         const double vh = vbar * h;
         const double term = std::exp(-odnf * vh);
         vscal1 *= term;
         eksum1 *= (term * term);
         lensc1 *= std::exp(vh);
         vbox1 *= std::exp(D * vh);


         // BBK vbar 2/2
         DelP = odnf * 2 * eksum1 - tr_vir;
         DelP = DelP - D * vbox1 * bath::atmsph / units::prescon;
         DelP_m = DelP / qbar;
         vbar = (DelP_m * h_2 + vbar + sd * R) / opgh2;


         // save
         eksum0 = eksum1;
         vbox0 = vbox1;
         vscal0 = vscal1;
         lensc0 = lensc1;
      }
   }


   // scale atomic velocities
   darray::scale(g::q0, n, vscal0, vx);
   darray::scale(g::q0, n, vscal0, vy);
   darray::scale(g::q0, n, vscal0, vz);
   double vscal2 = vscal0 * vscal0;
   eksum *= vscal2;
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         ekin[i][j] *= vscal2;


   // update the periodic box size
   lvec1 *= lensc0;
   lvec2 *= lensc0;
   lvec3 *= lensc0;
   set_default_recip_box();
}
}


void vv_lpiston_npt(int istep, time_prec dt)
{
   int vers1 = rc_flag & calc::vmask;
   bool save = !(istep % inform::iwrite);
   if (!save)
      vers1 &= ~calc::energy;


   const double R = normal<double>(); // N(0,1)
   const time_prec dt_2 = 0.5 * dt;


   // iL1 1/2
   lp(dt, R);


   // iLf 1/2
   propagate_velocity(dt_2, gx, gy, gz);


   // iLr
   const double term = vbar * dt_2;
   const double term2 = term * term;
   const double expterm = std::exp(term);
   const double eterm2 = expterm * expterm;
   constexpr double e2 = 1.0 / 6;
   constexpr double e4 = 1.0 / 120;
   constexpr double e6 = 1.0 / 5040;
   constexpr double e8 = 1.0 / 362880;
   // sinh(x)/x: Taylor series
   double poly = 1 + term2 * (e2 + term2 * (e4 + term2 * (e6 + term2 * e8)));
   poly *= expterm * dt;
   propagate_pos_axbv(eterm2, poly);
   copy_pos_to_xyz(true);


   energy(vers1);


   // iLf 2/2
   propagate_velocity(dt_2, gx, gy, gz);


   // iL1 2/2
   lp(dt, R);
}
}
