#include "energy.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdintg.h"
#include "mdpq.h"
#include "mdpt.h"
#include "random.h"
#include "rattle.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>


namespace tinker {
void velocity_verlet(int istep, time_prec dt_ps)
{
   int vers0 = rc_flag & calc::vmask;
   int vers1 = vers0;

   bool save = !(istep % inform::iwrite);
   bool mcbaro = false;
   if (barostat == MONTE_CARLO_BAROSTAT) {
      // toggle off the calc::virial bit if Monte Carlo Barostat is in use
      vers1 &= ~calc::virial;
      double rdm = random<double>();
      if (rdm < 1.0 / bath::voltrial)
         mcbaro = true;
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
      darray::copy(PROCEED_NEW_Q, n, rattle_xold, xpos);
      darray::copy(PROCEED_NEW_Q, n, rattle_yold, ypos);
      darray::copy(PROCEED_NEW_Q, n, rattle_zold, zpos);
   }
   // s += v * dt
   propagate_pos(dt_ps);
   if (userat)
      rattle(dt_ps, rattle_xold, rattle_yold, rattle_zold);
   copy_pos_to_xyz(true);

   // update gradient
   energy(vers1);

   // half-step corrections for certain thermostats and barostats
   halftime_correction(mcbaro);

   // v += a * dt/2
   propagate_velocity(dt_2, gx, gy, gz);
   if (userat)
      rattle2(dt_ps, vers1 bitand calc::virial);

   // full-step corrections
   T_prec temp;
   temper(dt_ps, temp);
}
}
