#include "energy.h"
#include "mathfunc_pow2.h"
#include "mdcalc.h"
#include "mddebug.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include "mdpt.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>


namespace tinker {
grad_prec *gx1, *gy1, *gz1;
grad_prec *gx2, *gy2, *gz2;


const TimeScaleConfig& respa_tsconfig()
{
   constexpr int fast = floor_log2_constexpr(RESPA_FAST); // short-range
   constexpr int slow = floor_log2_constexpr(RESPA_SLOW); // long-range
   static TimeScaleConfig tsconfig{
      {"ebond", fast},         {"eangle", fast},        {"estrbnd", fast},
      {"eurey", fast},         {"eopbend", fast},       {"etors", fast},
      {"eimprop", fast},       {"eimptor", fast},       {"epitors", fast},
      {"estrtor", fast},       {"eangtor", fast},       {"etortor", fast},
      {"egeom", fast},

      {"evalence", fast},

      {"evdw", slow},

      {"echarge", slow},       {"echglj", slow},

      {"emplar", slow},        {"empole", slow},        {"epolar", slow},

      {"empole_chgpen", slow}, {"epolar_chgpen", slow},

      {"echgtrn", slow},       {"edisp", slow},         {"erepel", slow},
      {"ehippo", slow},   
			
			{"empole_chgpen_aplus", slow},
			{"echgtrn_aplus", slow}, {"eaplus", slow}, 
   };
   return tsconfig;
}


/**
 * dT = nrsp dt
 * nrsp >= 1
 *
 * v += a_slow dT/2
 * DO irsp = 1, nrsp
 *    v += a_fast dt/2
 *    s += v dt
 *    update a_fast
 *    v += a_fast dt/2
 * END DO
 * update a_slow
 * v += a_slow dT/2
 * thermostat
 *
 * e.g. nrsp = 3
 * [v += a_slow dT/2]
 *    [v += a_fast dt/2] [s += v dt] [update a_fast] [v += a_fast dt/2]
 *    [v += a_fast dt/2] [s += v dt] [update a_fast] [v += a_fast dt/2]
 *    [v += a_fast dt/2] [s += v dt] [update a_fast] [v += a_fast dt/2]
 * [update a_slow] [v += a_slow dT/2]
 *
 * is equivalent to
 * [v += a_slow dT/2] [v += a_fast dt/2]
 *    [s += v dt] [update a_fast] [v += a_fast dt/2] [v += a_fast dt/2]
 *    [s += v dt] [update a_fast] [v += a_fast dt/2] [v += a_fast dt/2]
 * [s += v dt] [update a_fast] [v += a_fast dt/2]
 * [update a_slow] [v += a_slow dT/2]
 *
 * that is
 * [v += a_slow dT/2] [v += a_fast dt/2] ==> [v += (a_slow dT/2 + a_fast dt/2)]
 *    [s += v dt] [update a_fast] [v += a_fast dt]
 *    [s += v dt] [update a_fast] [v += a_fast dt]
 * [s += v dt]
 * [update a_fast] [update a_slow]
 * [v += a_fast dt/2] [v += a_slow dT/2] ==> [v += (a_fast dt/2 + a_slow dT/2)]
 */
void respa_fast_slow(int istep, time_prec dt_ps)
{
   int vers0 = rc_flag & calc::vmask;
   int vers1 = vers0;

   bool save = !(istep % inform::iwrite);
   bool mcbaro = do_pmonte;
   if (barostat == MONTE_CARLO_BAROSTAT) {
      // toggle off the calc::virial bit if Monte Carlo Barostat is in use
      vers1 &= ~calc::virial;
   }
   // toggle off the calc::energy bit if neither save nor mcbaro
   if (!save && !mcbaro)
      vers1 &= ~calc::energy;


   time_prec arespa = mdstuf::arespa;   // inner time step
   const time_prec eps = 1.0 / 1048576; // 2**-20
   int nalt = (int)(dt_ps / (arespa + eps)) + 1;
   time_prec dt_2 = 0.5 * dt_ps;
   time_prec dta = dt_ps / nalt;
   time_prec dta_2 = 0.5 * dta;


   virial_prec vir_fast[9] = {0};
   virial_prec vir_f[9];
   energy_prec esum_f;


   // g1: fast gradients; g2: slow gradients; see integrate_data()
   // v += a_slow dT/2
   // propagate_velocity(dt_2, gx2, gy2, gz2);
   // v += a_fast dt/2
   // propagate_velocity(dta_2, gx1, gy1, gz1);
   propagate_velocity2(dta_2, gx1, gy1, gz1, dt_2, gx2, gy2, gz2);


   for (int ifast = 1; ifast < nalt; ++ifast) {
      // s += v dt
      propagate_pos(dta);
      copy_pos_to_xyz(false);


      // update a_fast
      energy(vers1, RESPA_FAST, respa_tsconfig());
      copy_virial(vers1, vir_f);
      if (vers1 & calc::virial) {
         for (int i = 0; i < 9; ++i)
            vir_fast[i] += vir_f[i];
      }


      // v += a_fast dt
      propagate_velocity(dta, gx, gy, gz);
   }


   // s += v dt
   propagate_pos(dta);
   copy_pos_to_xyz(true);


   // update a_fast
   energy(vers1, RESPA_FAST, respa_tsconfig());
   darray::copy(g::q0, n, gx1, gx);
   darray::copy(g::q0, n, gy1, gy);
   darray::copy(g::q0, n, gz1, gz);
   copy_energy(vers1, &esum_f);
   copy_virial(vers1, vir_f);
   if (vers1 & calc::virial) {
      for (int i = 0; i < 9; ++i)
         vir_fast[i] += vir_f[i];
   }


   // update a_slow
   energy(vers1, RESPA_SLOW, respa_tsconfig());
   darray::copy(g::q0, n, gx2, gx);
   darray::copy(g::q0, n, gy2, gy);
   darray::copy(g::q0, n, gz2, gz);
   // esum: e slow
   // vir: v slow
   // esum_f: e fast
   // vir_fast: nalt total v fast
   if (vers1 & calc::energy)
      esum += esum_f;
   if (vers1 & calc::virial) {
      for (int i = 0; i < 9; ++i)
         vir[i] += vir_fast[i] / nalt;
   }


   // half-step corrections for certain thermostats and barostats
   T_prec temp;
   temper2(dt_ps, temp);
   pressure2(esum, temp);


   // v += a_fast dt/2
   // propagate_velocity(dta_2, gx1, gy1, gz1);
   // v += a_slow dT/2
   // propagate_velocity(dt_2, gx2, gy2, gz2);
   propagate_velocity2(dta_2, gx1, gy1, gz1, dt_2, gx2, gy2, gz2);


   // full-step corrections
   temper(dt_ps, temp, save);
   pressure(dt_ps);
}
}
