#include "mdpt.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "platform.h"
#include <cassert>


namespace tinker {
void kinetic(T_prec& temp)
{
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      kinetic_cu(temp);
   else
#endif
      kinetic_acc(temp);
}


void temper(time_prec dt, T_prec& temp)
{
   kinetic(temp);
   if (thermostat == NONE_THERMOSTAT)
      return;

   if (thermostat == BUSSI_THERMOSTAT)
      bussi_thermostat(dt, temp);
   else
      assert(false);
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


//====================================================================//


void pressure2(energy_prec epot, T_prec temp)
{
   if (do_pmonte) {
      monte_carlo_barostat(epot, temp);
   }
}
}
