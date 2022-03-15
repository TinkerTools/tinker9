#include "box.h"
#include "energy.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdpq.h"
#include "mdpt.h"
#include "rattle.h"
#include <tinker/detail/inform.hh>

namespace tinker {
void velocity_verlet(int istep, time_prec dt_ps)
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

   time_prec dt_2 = 0.5f * dt_ps;

   // gradients were calculated in integrate_data()
   // v += a * dt/2
   propagate_velocity(dt_2, gx, gy, gz);

   const bool userat = use_rattle();
   if (userat) {
      darray::copy(g::q0, n, rattle_xold, xpos);
      darray::copy(g::q0, n, rattle_yold, ypos);
      darray::copy(g::q0, n, rattle_zold, zpos);
   }
   // s += v * dt
   propagate_pos(dt_ps);
   if (userat)
      rattle(dt_ps, rattle_xold, rattle_yold, rattle_zold);
   copy_pos_to_xyz(true);

   // update gradient
   energy(vers1);

   // half-step corrections for certain thermostats and barostats
   T_prec temp;
   temper2(dt_ps, temp);
   pressure2(esum, temp);

   // v += a * dt/2
   propagate_velocity(dt_2, gx, gy, gz);
   if (userat)
      rattle2(dt_ps, vers1 & calc::virial);

   // full-step corrections
   temper(dt_ps, temp, save);
   pressure(dt_ps);
}
}
