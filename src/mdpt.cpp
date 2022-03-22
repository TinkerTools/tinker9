#include "md/md.h"
#include "platform.h"
#include <cassert>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>

namespace tinker {
void mdKinetic(T_prec& temp)
{
   mdKineticExplicit(temp, eksum, ekin, vx, vy, vz);
}

extern void kinetic_energy_acc(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n,
   const double* mass, const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);
extern void kinetic_energy_cu(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n,
   const double* mass, const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);
void mdKineticEnergy(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n,
   const double* mass, const vel_prec* vx, const vel_prec* vy, const vel_prec* vz)
{
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      kinetic_energy_cu(eksum_out, ekin_out, n, mass, vx, vy, vz);
   else
#endif
      kinetic_energy_acc(eksum_out, ekin_out, n, mass, vx, vy, vz);
}

void mdKineticExplicit(T_prec& temp_out, energy_prec& eksum_out, energy_prec (&ekin_out)[3][3],
   const vel_prec* vx, const vel_prec* vy, const vel_prec* vz)
{
   mdKineticEnergy(eksum_out, ekin_out, n, mass, vx, vy, vz);
   temp_out = 2 * eksum_out / (mdstuf::nfree * units::gasconst);
}

void mdBussiThermostat(time_prec dt, T_prec temp)
{
   mdBussiThermostat_acc(dt, temp);
}

//====================================================================//

void mdMonteCarloBarostat(energy_prec epot, T_prec temp)
{
   mdMonteCarloBarostat_acc(epot, temp);
}

void mdBerendsenBarostat(time_prec dt)
{
   mdBerendsenBarostat_acc(dt);
}
}
