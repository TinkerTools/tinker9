#include "md/pt.h"
#include "md/pq.h"
#include "tool/externfunc.h"
#include "tool/platform.h"
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>

namespace tinker {
TINKER_F2VOID(cu, 1, acc, 1, kineticEnergy, energy_prec&, energy_prec (&)[3][3], int n,
   const double*, const vel_prec*, const vel_prec*, const vel_prec*);
void kineticEnergy(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n, const double* mass,
   const vel_prec* vx, const vel_prec* vy, const vel_prec* vz)
{
   TINKER_F2CALL(cu, 1, acc, 1, kineticEnergy, eksum_out, ekin_out, n, mass, vx, vy, vz);
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
TINKER_F2VOID(cu, 0, acc, 1, bussiThermostat, time_prec, T_prec);
void bussiThermostat(time_prec dt, T_prec temp)
{
   TINKER_F2CALL(cu, 0, acc, 1, bussiThermostat, dt, temp);
}

TINKER_F2VOID(cu, 0, acc, 1, monteCarloBarostat, energy_prec, T_prec);
void monteCarloBarostat(energy_prec epot, T_prec temp)
{
   TINKER_F2CALL(cu, 0, acc, 1, monteCarloBarostat, epot, temp);
}

TINKER_F2VOID(cu, 0, acc, 1, berendsenBarostat, time_prec);
void berendsenBarostat(time_prec dt)
{
   TINKER_F2CALL(cu, 0, acc, 1, berendsenBarostat, dt);
}
}
