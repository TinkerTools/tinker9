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


void kinetic_leapfrog(T_prec& temp)
{
   // Ek at +1/2
   T_prec t1;
   energy_prec ekin1[3][3];
   kinetic_explicit(t1, eksum, ekin1, leapfrog_vx, leapfrog_vy, leapfrog_vz);


   // Ek at -1/2
   T_prec t2;
   energy_prec ekin2[3][3];
   kinetic_explicit(t2, eksum_old, ekin2, leapfrog_vxold, leapfrog_vyold,
                    leapfrog_vzold);


   ekin[0][0] = 0.5 * (ekin1[0][0] + ekin2[0][0]);
   ekin[0][1] = 0.5 * (ekin1[0][1] + ekin2[0][1]);
   ekin[0][2] = 0.5 * (ekin1[0][2] + ekin2[0][2]);
   ekin[1][0] = 0.5 * (ekin1[1][0] + ekin2[1][0]);
   ekin[1][1] = 0.5 * (ekin1[1][1] + ekin2[1][1]);
   ekin[1][2] = 0.5 * (ekin1[1][2] + ekin2[1][2]);
   ekin[2][0] = 0.5 * (ekin1[2][0] + ekin2[2][0]);
   ekin[2][1] = 0.5 * (ekin1[2][1] + ekin2[2][1]);
   ekin[2][2] = 0.5 * (ekin1[2][2] + ekin2[2][2]);
   temp = 0.5 * (t1 + t2);
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
