#include "md_pt.h"
#include "md_calc.h"
#include "md_egv.h"
#include "platform.h"
#include <cassert>


TINKER_NAMESPACE_BEGIN
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


void monte_carlo_barostat(energy_prec epot)
{
   monte_carlo_barostat_acc(epot);
}


//====================================================================//


void halftime_correction(bool do_voltrial)
{
   if (thermostat == NOSE_HOOVER_CHAIN_THERMOSTAT &&
       barostat == MONTE_CARLO_BAROSTAT) {
   } else if (thermostat == NOSE_HOOVER_CHAIN_THERMOSTAT) {
   } else if (barostat == MONTE_CARLO_BAROSTAT && do_voltrial) {
      energy_prec epot;
      copy_energy(calc::energy, &epot);
      monte_carlo_barostat(epot);
   }
}
TINKER_NAMESPACE_END
