#include "mdpt.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include "mdprec.h"
#include "platform.h"
#include <cassert>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>

namespace tinker {
void kinetic(T_prec& temp)
{
   kinetic_explicit(temp, eksum, ekin, vx, vy, vz);
}

extern void kinetic_energy_acc(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n,
   const double* mass, const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);
extern void kinetic_energy_cu(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n,
   const double* mass, const vel_prec* vx, const vel_prec* vy, const vel_prec* vz);
void kinetic_energy(energy_prec& eksum_out, energy_prec (&ekin_out)[3][3], int n,
   const double* mass, const vel_prec* vx, const vel_prec* vy, const vel_prec* vz)
{
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      kinetic_energy_cu(eksum_out, ekin_out, n, mass, vx, vy, vz);
   else
#endif
      kinetic_energy_acc(eksum_out, ekin_out, n, mass, vx, vy, vz);
}

void kinetic_explicit(T_prec& temp_out, energy_prec& eksum_out, energy_prec (&ekin_out)[3][3],
   const vel_prec* vx, const vel_prec* vy, const vel_prec* vz)
{
   kinetic_energy(eksum_out, ekin_out, n, mass, vx, vy, vz);
   temp_out = 2 * eksum_out / (mdstuf::nfree * units::gasconst);
}

void bussi_thermostat(time_prec dt, T_prec temp)
{
   bussi_thermostat_acc(dt, temp);
}

//====================================================================//

pos_prec *x_pmonte, *y_pmonte, *z_pmonte;

void monte_carlo_barostat(energy_prec epot, T_prec temp)
{
   monte_carlo_barostat_acc(epot, temp);
}

void berendsen_barostat(time_prec dt)
{
   berendsen_barostat_acc(dt);
}
}
