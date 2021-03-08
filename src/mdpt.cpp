#include "mdpt.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include "platform.h"
#include <cassert>


namespace tinker {
void kinetic(T_prec& temp)
{
   kinetic_explicit(temp, eksum, ekin, vx, vy, vz);
}


void kinetic_explicit(T_prec& temp_out, energy_prec& eksum_out,
                      energy_prec (&ekin_out)[3][3], const vel_prec* vx,
                      const vel_prec* vy, const vel_prec* vz)
{
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      kinetic_explicit_cu(temp_out, eksum_out, ekin_out, vx, vy, vz);
   else
#endif
      kinetic_explicit_acc(temp_out, eksum_out, ekin_out, vx, vy, vz);
}


void temper(time_prec dt, T_prec& temp, bool save)
{
   if (thermostat == BUSSI_THERMOSTAT) {
      kinetic(temp);
      bussi_thermostat(dt, temp);
   } else {
      // We don't need temperature for NVE but we still compute it anyway.
      if ((thermostat != NONE_THERMOSTAT and barostat != NONE_BAROSTAT) or
          save) {
         kinetic(temp);
      }
   }
}


void pressure(time_prec dt)
{
   if (barostat == NONE_BAROSTAT)
      return;

   if (barostat == BERENDSEN_BAROSTAT)
      berendsen_barostat(dt);
}


Thermostat thermostat;


void bussi_thermostat(time_prec dt, T_prec temp)
{
   bussi_thermostat_acc(dt, temp);
}


//====================================================================//


Barostat barostat;


pos_prec *x_pmonte, *y_pmonte, *z_pmonte;
vel_prec *vx_pmonte, *vy_pmonte, *vz_pmonte;
bool do_pmonte;


void monte_carlo_barostat(energy_prec epot, T_prec temp)
{
   monte_carlo_barostat_acc(epot, temp);
}


void berendsen_barostat(time_prec dt)
{
   berendsen_barostat_acc(dt);
}


//====================================================================//


void pressure2(energy_prec epot, T_prec temp)
{
   if (do_pmonte) {
      monte_carlo_barostat(epot, temp);
   }
}
}
