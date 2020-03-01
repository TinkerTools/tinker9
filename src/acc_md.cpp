#include "io_print.h"
#include "md.h"
#include "tinker_rt.h"
#include "wait_queue.h"
#include <tinker/detail/bound.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>

TINKER_NAMESPACE_BEGIN
void kinetic_acc(T_prec& temp)
{
   const energy_prec ekcal_inv = 1.0 / units::ekcal;
   energy_prec exx = 0;
   energy_prec eyy = 0;
   energy_prec ezz = 0;
   energy_prec exy = 0;
   energy_prec eyz = 0;
   energy_prec ezx = 0;
   #pragma acc parallel loop independent async\
               deviceptr(mass,vx,vy,vz)
   for (int i = 0; i < n; ++i) {
      energy_prec term = 0.5f * mass[i] * ekcal_inv;
      exx += term * vx[i] * vx[i];
      eyy += term * vy[i] * vy[i];
      ezz += term * vz[i] * vz[i];
      exy += term * vx[i] * vy[i];
      eyz += term * vy[i] * vz[i];
      ezx += term * vz[i] * vx[i];
   }
   ekin[0][0] = exx;
   ekin[0][1] = exy;
   ekin[0][2] = ezx;
   ekin[1][0] = exy;
   ekin[1][1] = eyy;
   ekin[1][2] = eyz;
   ekin[2][0] = ezx;
   ekin[2][1] = eyz;
   ekin[2][2] = ezz;
   eksum = exx + eyy + ezz;
   temp = 2 * eksum / (mdstuf::nfree * units::gasconst);
}

void mdrest_acc(int istep)
{
   if (!mdstuf::dorest)
      return;
   if ((istep % mdstuf::irest) != 0)
      return;

   const energy_prec ekcal = units::ekcal;

   // zero out the total mass and overall linear velocity

   mass_prec totmass = 0;
   vel_prec vtot1 = 0, vtot2 = 0, vtot3 = 0;

   // compute linear velocity of the system center of mass

   #pragma acc parallel loop independent async\
               deviceptr(mass,vx,vy,vz)
   for (int i = 0; i < n; ++i) {
      mass_prec weigh = mass[i];
      totmass += weigh;
      vtot1 += vx[i] * weigh;
      vtot2 += vy[i] * weigh;
      vtot3 += vz[i] * weigh;
   }

   vtot1 /= totmass;
   vtot2 /= totmass;
   vtot3 /= totmass;

   // compute translational kinetic energy of overall system

   energy_prec etrans = vtot1 * vtot1 + vtot2 * vtot2 + vtot3 * vtot3;
   etrans *= 0.5f * totmass / ekcal;

   energy_prec erot = 0;
   pos_prec xtot = 0, ytot = 0, ztot = 0;
   vel_prec vang[3] = {0, 0, 0}; // angular momentum
   if (!bound::use_bounds) {

      // find the center of mass coordinates of the overall system
      // compute the angular momentum of the overall system

      vel_prec mang1 = 0, mang2 = 0, mang3 = 0;

      #pragma acc parallel loop independent async\
                  deviceptr(mass,xpos,ypos,zpos,vx,vy,vz)
      for (int i = 0; i < n; ++i) {
         mass_prec weigh = mass[i];
         xtot += xpos[i] * weigh;
         ytot += ypos[i] * weigh;
         ztot += zpos[i] * weigh;
         mang1 += (ypos[i] * vz[i] - zpos[i] * vy[i]) * weigh;
         mang2 += (zpos[i] * vx[i] - xpos[i] * vz[i]) * weigh;
         mang3 += (xpos[i] * vy[i] - ypos[i] * vx[i]) * weigh;
      }
      xtot /= totmass;
      ytot /= totmass;
      ztot /= totmass;
      mang1 -= (ytot * vtot3 - ztot * vtot2) * totmass;
      mang2 -= (ztot * vtot1 - xtot * vtot3) * totmass;
      mang3 -= (xtot * vtot2 - ytot * vtot1) * totmass;

      // calculate the moment of inertia tensor

      pos_prec xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;

      #pragma acc parallel loop independent async\
                  deviceptr(mass,xpos,ypos,zpos)
      for (int i = 0; i < n; ++i) {
         mass_prec weigh = mass[i];
         pos_prec xdel = xpos[i] - xtot;
         pos_prec ydel = ypos[i] - ytot;
         pos_prec zdel = zpos[i] - ztot;
         xx += xdel * xdel * weigh;
         xy += xdel * ydel * weigh;
         xz += xdel * zdel * weigh;
         yy += ydel * ydel * weigh;
         yz += ydel * zdel * weigh;
         zz += zdel * zdel * weigh;
      }

      double tensor[3][3];
      double eps = (n <= 2 ? 0.000001 : 0);
      tensor[0][0] = yy + zz + eps;
      tensor[0][1] = -xy;
      tensor[0][2] = -xz;
      tensor[1][0] = -xy;
      tensor[1][1] = xx + zz + eps;
      tensor[1][2] = -yz;
      tensor[2][0] = -xz;
      tensor[2][1] = -yz;
      tensor[2][2] = xx + yy + eps;

      int ndim = 3;
      TINKER_RT(invert)(&ndim, &tensor[0][0]);

      // compute angular velocity and rotational kinetic energy

      for (int i = 0; i < 3; ++i) {
         vang[i] =
            tensor[0][i] * mang1 + tensor[1][i] * mang2 + tensor[2][i] * mang3;
      }
      erot = vang[0] * mang1 + vang[1] * mang2 + vang[2] * mang3;
      erot *= (0.5f / ekcal);
   }

   // eliminate any translation of the overall system

   #pragma acc parallel loop independent async deviceptr(vx,vy,vz)
   for (int i = 0; i < n; ++i) {
      vx[i] -= vtot1;
      vy[i] -= vtot2;
      vz[i] -= vtot3;
   }

   // print the translational velocity of the overall system

   if (inform::debug) {
      print(
         stdout,
         " System Linear Velocity :  {:12.2f}{:12.2f}{:12.2f}\n Translational "
         "Kinetic Energy :{:10s}{:12.4f} Kcal/mole\n",
         vtot1, vtot2, vtot3, "", etrans);
   }

   // eliminate any rotation about the system center of mass

   if (!bound::use_bounds) {
      #pragma acc parallel loop independent async\
                  deviceptr(xpos,ypos,zpos,vx,vy,vz)
      for (int i = 0; i < n; ++i) {
         pos_prec xdel = xpos[i] - xtot;
         pos_prec ydel = ypos[i] - ytot;
         pos_prec zdel = zpos[i] - ztot;
         vx[i] -= vang[1] * zdel + vang[2] * ydel;
         vy[i] -= vang[2] * xdel + vang[0] * zdel;
         vz[i] -= vang[0] * ydel + vang[1] * xdel;
      }
   }

   // print the angular velocity of the overall system

   if (inform::debug) {
      print(stdout,
            " System Angular Velocity : {:12.2f}{:12.2f}{:12.2f}\n Rotational "
            "Kinetic Energy :{:13s}{:12.4f} Kcal/mole\n",
            vang[0], vang[1], vang[2], "", erot);
   }
}

void copy_pos_to_xyz_acc()
{
   if (sizeof(pos_prec) == sizeof(real))
      return;


   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,xpos,ypos,zpos)
   for (int i = 0; i < n; ++i) {
      x[i] = xpos[i];
      y[i] = ypos[i];
      z[i] = zpos[i];
   }
}

void propagate_pos_acc(time_prec dt)
{
   #pragma acc parallel loop independent async\
               deviceptr(xpos,ypos,zpos,vx,vy,vz)
   for (int i = 0; i < n; ++i) {
      xpos[i] += dt * vx[i];
      ypos[i] += dt * vy[i];
      zpos[i] += dt * vz[i];
   }
}

void propagate_velocity_acc(time_prec dt, const real* grx, const real* gry,
                            const real* grz)
{
   const vel_prec ekcal = units::ekcal;
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vx,vy,vz,grx,gry,grz)
   for (int i = 0; i < n; ++i) {
      vel_prec coef = -ekcal * massinv[i] * dt;
      vx[i] += coef * grx[i];
      vy[i] += coef * gry[i];
      vz[i] += coef * grz[i];
   }
}

void propagate_velocity2_acc(time_prec dt, const real* grx, const real* gry,
                             const real* grz, time_prec dt2, const real* grx2,
                             const real* gry2, const real* grz2)
{
   const vel_prec ekcal = units::ekcal;
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vx,vy,vz,grx,gry,grz,grx2,gry2,grz2)
   for (int i = 0; i < n; ++i) {
      vel_prec coef = -ekcal * massinv[i];
      vx[i] += coef * (grx[i] * dt + grx2[i] * dt2);
      vy[i] += coef * (gry[i] * dt + gry2[i] * dt2);
      vz[i] += coef * (grz[i] * dt + grz2[i] * dt2);
   }
}

void zero_gradient(DMFlag flag, size_t nelem, real* gx, real* gy, real* gz)
{
   bool sync = flag & DMFlag::DEFAULT_Q;
   if (sync) {
      #pragma acc parallel loop independent deviceptr(gx,gy,gz)
      for (int i = 0; i < nelem; ++i) {
         gx[i] = 0;
         gy[i] = 0;
         gz[i] = 0;
      }
   } else {
      #pragma acc parallel loop independent async deviceptr(gx,gy,gz)
      for (int i = 0; i < nelem; ++i) {
         gx[i] = 0;
         gy[i] = 0;
         gz[i] = 0;
      }
   }
   // if (flag & DMFlag::WAIT) {
   wait_queue(flag);
   // }
}
TINKER_NAMESPACE_END
