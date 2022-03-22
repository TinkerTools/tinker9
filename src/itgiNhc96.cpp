#include "ff/box.h"
#include "ff/energy.h"
#include "integrator.h"
#include "md.h"
#include "tool/error.h"
#include <cmath>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>

namespace tinker {
static double gbar;
static double vnh[maxnose];
static double qnh[maxnose];
static double gnh[maxnose];
static double press;

static void hoover(time_prec dt, virial_prec press)
{
   constexpr int nc = 5;
   constexpr int ns = 3;
   // w[0] = 1/(2 - 2**(1/3))
   // w[1] = 1 - w[0] - w[2]
   // w[2] = w[0]
   constexpr double w[3] = {1.351207191959657634047687808971460826921999376217144828328,
      -1.70241438391931526809537561794292165384399875243428965665,
      1.351207191959657634047687808971460826921999376217144828328};

   const double vbox = boxVolume();
   T_prec temp;
   mdKinetic(temp);
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
   darray::scale(g::q0, n, scale, vx);
   darray::scale(g::q0, n, scale, vy);
   darray::scale(g::q0, n, scale, vz);
}

static void nhc_npt(int istep, time_prec dt)
{
   int vers1 = rc_flag & calc::vmask;
   bool save = !(istep % inform::iwrite);
   if (!save)
      vers1 &= ~calc::energy;

   // set some time values for the dynamics integration
   const time_prec dt_2 = 0.5f * dt;
   // This initialization is intentially kept the same as the Fortran code.
   // Yes, the real press is available here, but when the Fortran code was
   // written, virial may not be available if simulation was restarted from
   // a ".dyn" file.
   if (istep == 1)
      press = bath::atmsph;

   // update thermostat and barostat values, scale atomic velocities
   hoover(dt, press);

   mdVel(dt_2, gx, gy, gz);

   double term = vbar * dt_2;
   double term2 = term * term;
   double expterm = std::exp(term);
   double eterm2 = expterm * expterm;

   // update the periodic box size and total volume
   // eq. 42 volume
   lvec1 *= eterm2;
   lvec2 *= eterm2;
   lvec3 *= eterm2;
   boxSetDefaultRecip();

   // update atomic positions via coupling to barostat
   // eq. 42 coordinates
   constexpr double e2 = 1.0 / 6;
   constexpr double e4 = e2 / 20;
   constexpr double e6 = e4 / 42;
   constexpr double e8 = e6 / 72;
   // sinh(x)/x: Taylor series
   double poly = 1 + term2 * (e2 + term2 * (e4 + term2 * (e6 + term2 * e8)));
   poly *= expterm * dt;
   mdPosAxbv(eterm2, poly);
   mdCopyPosToXyz(true);

   energy(vers1);

   mdVel(dt_2, gx, gy, gz);

   // update thermostat and barostat values, scale atomic velocities
   hoover(dt, press);

   // set isotropic pressure to the average of tensor diagonal
   double vbox = boxVolume();
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

namespace tinker {
const char* Nhc96Integrator::name() const
{
   return "Molecular Dynamics Trajectory via Nose-Hoover NPT Algorithm";
}

void Nhc96Integrator::kickoff()
{
   double ekt = units::gasconst * bath::kelvin;
   vbar = 0;
   qbar = (mdstuf::nfree + 1) * ekt * bath::taupres * bath::taupres;
   gbar = 0;
   for (int i = 0; i < maxnose; ++i) {
      vnh[i] = 0;
      qnh[i] = ekt * bath::tautemp * bath::tautemp;
      gnh[i] = 0;
   }
   qnh[0] *= mdstuf::nfree;
   energy(calc::grad | calc::virial);
}

Nhc96Integrator::Nhc96Integrator()
   : BasicIntegrator()
{
   if (useRattle())
      TINKER_THROW("Constraints under NH-NPT require the ROLL algorithm.");

   this->kickoff();
}

void Nhc96Integrator::dynamic(int istep, time_prec dt)
{
   nhc_npt(istep, dt);
}
}
