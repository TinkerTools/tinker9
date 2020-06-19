#include "nose.h"
#include "box.h"
#include "energy.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdpq.h"
#include "mdpt.h"
#include <cmath>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>


namespace tinker {
double vbar;
double qbar;
double gbar;
double vnh[maxnose];
double qnh[maxnose];
double gnh[maxnose];


void hoover(time_prec dt, virial_prec press)
{
   constexpr int nc = 5;
   constexpr int ns = 3;
   // w[0] = 1/(2 - 2**(1/3))
   // w[1] = 1 - w[0] - w[2]
   // w[2] = w[0]
   constexpr double w[3] = {
      1.351207191959657634047687808971460826921999376217144828328,
      -1.70241438391931526809537561794292165384399875243428965665,
      1.351207191959657634047687808971460826921999376217144828328};


   const double vbox = volbox();
   T_prec temp;
   kinetic(temp);
   const double ekt = units::gasconst * bath::kelvin;
   const int df = mdstuf::nfree;
   const double odnf = 1 + 3.0 / df;
   const double gn1kt = (1 + df) * ekt;
   const double dpress = (press - bath::atmsph) / units::prescon;
   const time_prec dtc = dt / nc;


   double scale = 1.0;
   for (int k = 0; k < nc; ++k) {
      for (int j = 0; j < ns; ++j) {
         const time_prec dts = w[j] * dtc;
         const time_prec dt2 = 0.5 * dts;
         const time_prec dt4 = 0.25 * dts;
         const time_prec dt8 = 0.125 * dts;
         double expterm;


         // update barostat and thermostat velocities and forces
         // eq. 41 first half
         for (int i = maxnose - 1; i > -1; --i) {
            if (i == 0)
               gnh[i] = (2 * eksum + qbar * vbar * vbar - gn1kt) / qnh[i];
            else
               gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - ekt) / qnh[i];


            if (i == maxnose - 1)
               vnh[i] += gnh[i] * dt4;
            else {
               double exptm = std::exp(-vnh[i + 1] * dt8);
               vnh[i] = (vnh[i] * exptm + gnh[i] * dt4) * exptm;
            }
         }
         gbar = (2 * eksum * odnf + 3 * vbox * dpress) / qbar;
         expterm = std::exp(-vnh[0] * dt8);
         vbar = (vbar * expterm + gbar * dt4) * expterm;


         // find velocity scale factor and update kinetic energy
         // eq. 41 velocities
         expterm = std::exp(-(vnh[0] + vbar * odnf) * dt2);
         scale *= expterm;
         double exptm2 = expterm * expterm;
         eksum *= exptm2;
         for (int ii = 0; ii < 3; ++ii)
            for (int jj = 0; jj < 3; ++jj)
               ekin[ii][jj] *= exptm2;


         // update barostat and thermostat velocities and forces
         // eq. 41 second half
         gbar = (2 * eksum * odnf + 3 * vbox * dpress) / qbar;
         expterm = std::exp(-vnh[0] * dt8);
         vbar = (vbar * expterm + gbar * dt4) * expterm;
         for (int i = 0; i < maxnose; ++i) {
            if (i == 0)
               gnh[i] = (2 * eksum + qbar * vbar * vbar - gn1kt) / qnh[i];
            else
               gnh[i] = (qnh[i - 1] * vnh[i - 1] * vnh[i - 1] - ekt) / qnh[i];


            if (i == maxnose - 1)
               vnh[i] += gnh[i] * dt4;
            else {
               double exptm = std::exp(-vnh[i + 1] * dt8);
               vnh[i] = (vnh[i] * exptm + gnh[i] * dt4) * exptm;
            }
         }
      }
   }


   // use scale factor to update the atomic velocities
   // eq. 41 velocities
   darray::scale(PROCEED_NEW_Q, n, scale, vx);
   darray::scale(PROCEED_NEW_Q, n, scale, vy);
   darray::scale(PROCEED_NEW_Q, n, scale, vz);
}


namespace {
double press;
}


void nhc_npt_xo_respa(int istep, time_prec dt)
{
   int vers1 = rc_flag & calc::vmask;
   bool save = !(!istep % inform::iwrite);
   if (!save)
      vers1 &= ~calc::energy;


   const time_prec arespa = mdstuf::arespa;  // inner time step
   constexpr time_prec eps = 1.0f / 1048576; // 2**-20
   const int nalt = (int)(dt / (arespa + eps)) + 1;
   const time_prec dt_2 = 0.5f * dt;
   const time_prec dta = dt / nalt;
   const time_prec dta_2 = 0.5f * dta;


   virial_prec vir_fast[9] = {0};
   virial_prec vir_f[9];
   energy_prec esum_f;


   // set some time values for the dynamics integration
   if (istep == 1)
      press = bath::atmsph;


   // update thermostat and barostat values, scale atomic velocities
   hoover(dt, press);


   propagate_velocity2(dta_2, gx1, gy1, gz1, dt_2, gx2, gy2, gz2);


   for (int ifast = 1; ifast <= nalt; ++ifast) {
      bool ilast = (ifast == nalt);


      double term = vbar * dta_2;
      double term2 = term * term;
      double expterm = std::exp(term);
      double eterm2 = expterm * expterm;


      // update the periodic box size and total volume
      // eq. 42 volume
      lvec1 *= eterm2;
      lvec2 *= eterm2;
      lvec3 *= eterm2;
      set_default_recip_box();


      // update atomic positions via coupling to barostat
      // eq. 42 coordinates
      constexpr double e2 = 1.0 / 6;
      constexpr double e4 = e2 / 20;
      constexpr double e6 = e4 / 42;
      constexpr double e8 = e6 / 72;
      // sinh(x)/x: Taylor series up to x**10
      double poly = 1 + term2 * (e2 + term2 * (e4 + term2 * (e6 + term2 * e8)));
      poly *= expterm * dta;
      propagate_xyz_axbv(eterm2, poly, ilast); // check nblist if ifast == nalt


      if (!ilast) {
         // update a_fast
         energy(vers1, RESPA_FAST, respa_tsconfig());
         copy_virial(vers1, vir_f);
         for (int i = 0; i < 9; ++i)
            vir_fast[i] += vir_f[i];


         // v += a_fast dt
         propagate_velocity(dta, gx, gy, gz);
      } else {
         // update a_fast
         energy(vers1, RESPA_FAST, respa_tsconfig());
         darray::copy(PROCEED_NEW_Q, n, gx1, gx);
         darray::copy(PROCEED_NEW_Q, n, gy1, gy);
         darray::copy(PROCEED_NEW_Q, n, gz1, gz);
         copy_energy(vers1, &esum_f);
         copy_virial(vers1, vir_f);
         for (int i = 0; i < 9; ++i)
            vir_fast[i] += vir_f[i];


         // update a_slow
         energy(vers1, RESPA_SLOW, respa_tsconfig());
         darray::copy(PROCEED_NEW_Q, n, gx2, gx);
         darray::copy(PROCEED_NEW_Q, n, gy2, gy);
         darray::copy(PROCEED_NEW_Q, n, gz2, gz);
         // esum: e slow
         // vir: v slow
         // esum_f: e fast
         // vir_fast: nalt total v fast
         if (vers1 & calc::energy)
            esum += esum_f;
         for (int i = 0; i < 9; ++i)
            vir[i] += vir_fast[i] / nalt;
      }
   }


   propagate_velocity2(dta_2, gx1, gy1, gz1, dt_2, gx2, gy2, gz2);


   // update thermostat and barostat values, scale atomic velocities
   hoover(dt, press);


   // set isotropic pressure to the average of tensor diagonal
   double vbox = volbox();
   double factor = units::prescon / vbox;
   double stress[3][3];
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         stress[i][j] = factor * (-vir[3 * i + j]);
      }
   }
   press = (stress[0][0] + stress[1][1] + stress[2][2]) / 3;
}
}
