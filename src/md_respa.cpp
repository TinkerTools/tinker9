#include "energy.h"
#include "mathfunc_pow2.h"
#include "md.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>


TINKER_NAMESPACE_BEGIN
const TimeScaleConfig& respa_tsconfig()
{
   constexpr int fast = floor_log2_constexpr(RESPA_FAST); // short-range
   constexpr int slow = floor_log2_constexpr(RESPA_SLOW); // long-range
   static TimeScaleConfig tsconfig{
      {"ebond", fast},   {"eangle", fast},    {"estrbnd", fast},
      {"eurey", fast},   {"eopbend", fast},   {"etors", fast},
      {"epitors", fast}, {"etortor", fast},   {"egeom", fast},

      {"evdw", slow},    {"elec_init", slow}, {"torque", slow},
      {"emplar", slow},  {"empole", slow},    {"epolar", slow},
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
void respa_fast_slow(int istep, real dt_ps)
{
   bool save = !(istep % inform::iwrite);
   int vers0 = rc_flag & (calc::virial | calc::grad | calc::energy);
   int vers1 = rc_flag & (calc::virial | calc::grad);


   real arespa = mdstuf::arespa;    // inner time step
   const real eps = 1.0f / 1048576; // 2**-20
   int nalt = (int)(dt_ps / (arespa + eps)) + 1;
   real dt_2 = 0.5f * dt_ps;
   real dta = dt_ps / nalt;
   real dta_2 = 0.5f * dta;


   real vir_fast[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
   real vir_f[9];
   real esum_f;


   // g1: fast gradients; g2: slow gradients; see integrate_data()
   // v += a_slow dT/2
   // propagate_velocity(dt_2, gx2, gy2, gz2);
   // v += a_fast dt/2
   // propagate_velocity(dta_2, gx1, gy1, gz1);
   propagate_velocity2(dta_2, gx1, gy1, gz1, dt_2, gx2, gy2, gz2);


   for (int ifast = 1; ifast < nalt; ++ifast) {
      // s += v dt
      propagate_xyz(dta, false);
      // update a_fast
      energy(vers1, RESPA_FAST, respa_tsconfig());
      copy_energy(vers1, nullptr, nullptr, nullptr, nullptr, vir_f);
      if (vers1 & calc::virial) {
         for (int i = 0; i < 9; ++i)
            vir_fast[i] += vir_f[i];
      }
      // v += a_fast dt
      propagate_velocity(dta, gx, gy, gz);
   }


   // s += v dt
   propagate_xyz(dta, true);
   // update a_fast
   if (save) {
      energy(vers0, RESPA_FAST, respa_tsconfig());
      copy_energy(vers0, &esum_f, gx1, gy1, gz1, vir_f);
   } else {
      energy(vers1, RESPA_FAST, respa_tsconfig());
      copy_energy(vers1, nullptr, gx1, gy1, gz1, vir_f);
   }
   if (rc_flag & calc::virial) {
      for (int i = 0; i < 9; ++i)
         vir_fast[i] += vir_f[i];
   }
   // update a_slow
   if (save) {
      energy(vers0, RESPA_SLOW, respa_tsconfig());
      copy_energy(vers0, nullptr, gx2, gy2, gz2, nullptr);
      if (rc_flag & calc::energy)
         esum += esum_f;
   } else {
      energy(vers1, RESPA_SLOW, respa_tsconfig());
      copy_energy(vers1, nullptr, gx2, gy2, gz2, nullptr);
   }
   if (rc_flag & calc::virial) {
      for (int i = 0; i < 9; ++i)
         vir[i] += vir_fast[i] / nalt;
   }
   // v += a_fast dt/2
   // propagate_velocity(dta_2, gx1, gy1, gz1);
   // v += a_slow dT/2
   // propagate_velocity(dt_2, gx2, gy2, gz2);
   propagate_velocity2(dta_2, gx1, gy1, gz1, dt_2, gx2, gy2, gz2);


   real temp;
   temper(dt_ps, temp);
}
TINKER_NAMESPACE_END
