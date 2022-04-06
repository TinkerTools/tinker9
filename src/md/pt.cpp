#include "md/pt.h"
#include "md/pq.h"
#include "tool/platform.h"
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>

namespace tinker {
extern void kineticEnergy_acc(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n,
   const double* mass, const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);
extern void kineticEnergy_cu(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n,
   const double* mass, const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);
void kineticEnergy(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n, const double* mass,
   const vel_prec* vx, const vel_prec* vy, const vel_prec* vz)
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      kineticEnergy_cu(eksum_out, ekin_out, n, mass, vx, vy, vz);
   else
#endif
      kineticEnergy_acc(eksum_out, ekin_out, n, mass, vx, vy, vz);
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
extern void bussiThermostat_acc(time_prec dt, T_prec temp);
void bussiThermostat(time_prec dt, T_prec temp)
{
   bussiThermostat_acc(dt, temp);
}

extern void monteCarloBarostat_acc(energy_prec epot, T_prec temp);
void monteCarloBarostat(energy_prec epot, T_prec temp)
{
   monteCarloBarostat_acc(epot, temp);
}

extern void berendsenBarostat_acc(time_prec);
void berendsenBarostat(time_prec dt)
{
   berendsenBarostat_acc(dt);
}
}
