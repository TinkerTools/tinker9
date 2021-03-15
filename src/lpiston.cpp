#include "lpiston.h"
#include "box.h"
#include "energy.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include "mdpt.h"
#include "nose.h"
#include "random.h"
#include "rattle.h"
#include <cmath>
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/stodyn.hh>
#include <tinker/detail/units.hh>


namespace tinker {
double lp_alpha;
double lp_molpres;
double* lp_molpres_buf;


void lp_molpressure(double alpha, double& val)
{
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      lp_molpressure_cu(alpha, val);
   else
#endif
      lp_molpressure_acc(alpha, val);
}


void lprat(time_prec dt, const pos_prec* xold, const pos_prec* yold,
           const pos_prec* zold)
{
   lprat_settle_acc(dt, xold, yold, zold);
   lprat_ch_acc(dt, xold, yold, zold);
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      lprat_methyl_cu(dt, xold, yold, zold);
#endif
   lprat_acc(dt, xold, yold, zold);
}


void lprat2(time_prec dt)
{
   lprat2_settle_acc(dt);
   lprat2_ch_acc(dt);
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      lprat2_methyl_cu(dt);
#endif
   lprat2_acc(dt);
}


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
   const double odnf = lp_alpha;


   const int nc = 2;
   const double h = dt / (2 * nc);
   const double h_2 = 0.5 * h;
   const double h_4 = 0.25 * h;
   const double opgh2 = 1.0 + g * h_2;
   const double omgh2 = 1.0 - g * h_2;
   const double sdbar = std::sqrt(2.0 * kbt * g * dt / qbar) / (4 * nc);


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
      for (int i = maxnose - 1; i > -1; --i) {
         if (i == 0)
            gnh[i] = (2 * eksum1 - Nf * kbt) / qnh[i];
         else
            gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - kbt) / qnh[i];


         if (i == maxnose - 1)
            vnh[i] += gnh[i] * h_2;
         else {
            double exptm = std::exp(-vnh[i + 1] * h_4);
            vnh[i] = (vnh[i] * exptm + gnh[i] * h_2) * exptm;
         }
      }


      // vbar 1/2
      DelP = (odnf * 2 * eksum1 - tr_vir);
      DelP = DelP - D * vol0 * bath::atmsph / units::prescon;
      gbar = DelP / qbar;
      vbar = omgh2 * vbar + gbar * h_2 + sdbar * R;


      // velocity
      double scal = std::exp(-h * (odnf * vbar + vnh[0]));
      velsc1 *= scal;
      eksum1 *= (scal * scal);


      // vbar 2/2
      DelP = (odnf * 2 * eksum1 - tr_vir);
      DelP = DelP - D * vol0 * bath::atmsph / units::prescon;
      gbar = DelP / qbar;
      vbar = (vbar + gbar * h_2 + sdbar * R) / opgh2;


      // vnh 2/2
      for (int i = 0; i < maxnose; ++i) {
         if (i == 0)
            gnh[i] = (2 * eksum1 - Nf * kbt) / qnh[i];
         else
            gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - kbt) / qnh[i];


         if (i == maxnose - 1)
            vnh[i] += gnh[i] * h_2;
         else {
            double exptm = std::exp(-vnh[i + 1] * h_4);
            vnh[i] = (vnh[i] * exptm + gnh[i] * h_2) * exptm;
         }
      }


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
}


void vv_lpiston_uc(int istep, time_prec dt)
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


void vv_lpiston_npt(int istep, time_prec dt)
{
   if (use_rattle()) {
      vv_lpiston_hc_acc(istep, dt);
   } else {
      vv_lpiston_uc(istep, dt);
   }
}
}
