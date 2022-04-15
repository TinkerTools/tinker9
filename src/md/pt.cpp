#include "math/random.h"
#include "md/misc.h"
#include "md/pq.h"
#include "tool/externfunc.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>

#include <cmath>

namespace tinker {
TINKER_FVOID2(cu, 1, acc, 1, kineticEnergy, energy_prec&, energy_prec (&)[3][3], int n,
   const double*, const vel_prec*, const vel_prec*, const vel_prec*);
void kineticEnergy(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n, const double* mass,
   const vel_prec* vx, const vel_prec* vy, const vel_prec* vz)
{
   TINKER_FCALL2(cu, 1, acc, 1, kineticEnergy, eksum_out, ekin_out, n, mass, vx, vy, vz);
}

void kineticExplicit(T_prec& temp_out, energy_prec& eksum_out, energy_prec (&ekin_out)[3][3],
   const vel_prec* vx, const vel_prec* vy, const vel_prec* vz)
{
   kineticEnergy(eksum_out, ekin_out, n, mass, vx, vy, vz);
   temp_out = 2 * eksum_out / (mdstuf::nfree * units::gasconst);
}

void kinetic(T_prec& temp)
{
   kineticExplicit(temp, eksum, ekin, vx, vy, vz);
}
}

namespace tinker {
void bussiThermostat(time_prec dt_prec, T_prec temp_prec)
{
   double dt = dt_prec;
   double temp = temp_prec;

   double tautemp = bath::tautemp;
   double kelvin = bath::kelvin;
   int nfree = mdstuf::nfree;
   double& eta = bath::eta;

   if (temp == 0)
      temp = 0.1;

   double c = std::exp(-dt / tautemp);
   double d = (1 - c) * (kelvin / temp) / nfree;
   double r = normal<double>();
   double s = chiSquared<double>(nfree - 1);
   double scale = c + (s + r * r) * d + 2 * r * std::sqrt(c * d);
   scale = std::sqrt(scale);
   if (r + std::sqrt(c / d) < 0)
      scale = -scale;
   eta *= scale;

   vel_prec sc = scale;
   mdVelScale(sc, n, vx, vy, vz);
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         ekin[i][j] *= scale * scale;
      }
   }
}

TINKER_FVOID2(cu, 0, acc, 1, monteCarloBarostat, energy_prec, T_prec);
void monteCarloBarostat(energy_prec epot, T_prec temp)
{
   TINKER_FCALL2(cu, 0, acc, 1, monteCarloBarostat, epot, temp);
}

TINKER_FVOID2(cu, 0, acc, 1, berendsenBarostat, time_prec);
void berendsenBarostat(time_prec dt)
{
   TINKER_FCALL2(cu, 0, acc, 1, berendsenBarostat, dt);
}
}
