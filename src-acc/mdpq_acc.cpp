#include "add.h"
#include "glob.molecule.h"
#include "image.h"
#include "mdpq.h"
#include <tinker/detail/units.hh>


namespace tinker {
void copy_pos_to_xyz_acc()
{
   if CONSTEXPR (sizeof(pos_prec) == sizeof(real))
      return;
   else
      #pragma acc parallel loop independent async\
                  deviceptr(x,y,z,xpos,ypos,zpos)
      for (int i = 0; i < n; ++i) {
         x[i] = xpos[i];
         y[i] = ypos[i];
         z[i] = zpos[i];
      }
}


void propagate_pos_acc(time_prec dt, pos_prec* qx, pos_prec* qy, pos_prec* qz,
                       const vel_prec* vlx, const vel_prec* vly,
                       const vel_prec* vlz)
{
   #pragma acc parallel loop independent async\
               deviceptr(qx,qy,qz,vlx,vly,vlz)
   for (int i = 0; i < n; ++i) {
      qx[i] += dt * vlx[i];
      qy[i] += dt * vly[i];
      qz[i] += dt * vlz[i];
   }
}


void propagate_pos_axbv_acc(double a, double b)
{
   pos_prec sa = a, sb = b;
   #pragma acc parallel loop independent async\
               deviceptr(xpos,ypos,zpos,vx,vy,vz)
   for (int i = 0; i < n; ++i) {
      xpos[i] = sa * xpos[i] + sb * vx[i];
      ypos[i] = sa * ypos[i] + sb * vy[i];
      zpos[i] = sa * zpos[i] + sb * vz[i];
   }
}


void bounds_pos_acc()
{
   auto nmol = molecule.nmol;
   const auto* imol = molecule.imol;
   const auto* kmol = molecule.kmol;
   const auto* molmass = molecule.molmass;


   #pragma acc parallel loop independent async\
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
               deviceptr(imol,kmol,xpos,ypos,zpos,mass,molmass)
   for (int i = 0; i < nmol; ++i) {
      // locate the center of each molecule
      pos_prec xmid = 0, ymid = 0, zmid = 0;
      int start = imol[i][0];
      int stop = imol[i][1];
      #pragma acc loop seq
      for (int j = start; j < stop; ++j) {
         int k = kmol[j];
         xmid += xpos[k] * mass[k];
         ymid += ypos[k] * mass[k];
         zmid += zpos[k] * mass[k];
      }
      auto weigh = molmass[i];
      xmid /= weigh;
      ymid /= weigh;
      zmid /= weigh;


      // locate the image of the center inside PBC box
      real xc, yc, zc;
      xc = xmid;
      yc = ymid;
      zc = zmid;
      image(xc, yc, zc);


      #pragma acc loop seq
      for (int j = start; j < stop; ++j) {
         int k = kmol[j];
         xpos[k] += xc - xmid;
         ypos[k] += yc - ymid;
         zpos[k] += zc - zmid;
      }
   }
}


//====================================================================//


void propagate_velocity_acc(time_prec dt, vel_prec* vlx, vel_prec* vly,
                            vel_prec* vlz, const vel_prec* vlx0,
                            const vel_prec* vly0, const vel_prec* vlz0,
                            const grad_prec* grx, const grad_prec* gry,
                            const grad_prec* grz)
{
   const vel_prec ekcal = units::ekcal;
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vlx,vly,vlz,vlx0,vly0,vlz0,grx,gry,grz)
   for (int i = 0; i < n; ++i) {
      vel_prec coef = -ekcal * massinv[i] * dt;
#if TINKER_DETERMINISTIC_FORCE
      vlx[i] = vlx0[i] + coef * to_flt_acc<vel_prec>(grx[i]);
      vly[i] = vly0[i] + coef * to_flt_acc<vel_prec>(gry[i]);
      vlz[i] = vlz0[i] + coef * to_flt_acc<vel_prec>(grz[i]);
#else
      vlx[i] = vlx0[i] + coef * grx[i];
      vly[i] = vly0[i] + coef * gry[i];
      vlz[i] = vlz0[i] + coef * grz[i];
#endif
   }
}


void propagate_velocity_acc(time_prec dt, vel_prec* vlx, vel_prec* vly,
                            vel_prec* vlz, const grad_prec* grx,
                            const grad_prec* gry, const grad_prec* grz)
{
   const vel_prec ekcal = units::ekcal;
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vlx,vly,vlz,grx,gry,grz)
   for (int i = 0; i < n; ++i) {
      vel_prec coef = -ekcal * massinv[i] * dt;
#if TINKER_DETERMINISTIC_FORCE
      vlx[i] += coef * to_flt_acc<vel_prec>(grx[i]);
      vly[i] += coef * to_flt_acc<vel_prec>(gry[i]);
      vlz[i] += coef * to_flt_acc<vel_prec>(grz[i]);
#else
      vlx[i] += coef * grx[i];
      vly[i] += coef * gry[i];
      vlz[i] += coef * grz[i];
#endif
   }
}


void propagate_velocity2_acc(time_prec dt, const grad_prec* grx,
                             const grad_prec* gry, const grad_prec* grz,
                             time_prec dt2, const grad_prec* grx2,
                             const grad_prec* gry2, const grad_prec* grz2)
{
   const vel_prec ekcal = units::ekcal;
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vx,vy,vz,grx,gry,grz,grx2,gry2,grz2)
   for (int i = 0; i < n; ++i) {
      vel_prec coef = -ekcal * massinv[i];
#if TINKER_DETERMINISTIC_FORCE
      // clang-format off
      vx[i] += coef*(to_flt_acc<vel_prec>(grx[i])*dt+to_flt_acc<vel_prec>(grx2[i])*dt2);
      vy[i] += coef*(to_flt_acc<vel_prec>(gry[i])*dt+to_flt_acc<vel_prec>(gry2[i])*dt2);
      vz[i] += coef*(to_flt_acc<vel_prec>(grz[i])*dt+to_flt_acc<vel_prec>(grz2[i])*dt2);
      // clang-format on
#else
      vx[i] += coef * (grx[i] * dt + grx2[i] * dt2);
      vy[i] += coef * (gry[i] * dt + gry2[i] * dt2);
      vz[i] += coef * (grz[i] * dt + grz2[i] * dt2);
#endif
   }
}


//====================================================================//


void swap_velocity_acc(vel_prec* vxnew, vel_prec* vynew, vel_prec* vznew,
                       vel_prec* vxold, vel_prec* vyold, vel_prec* vzold)
{
   #pragma acc parallel loop independent async\
               deviceptr(vxnew,vynew,vznew,vxold,vyold,vzold)
   for (int i = 0; i < n; ++i) {
      vel_prec sx, sy, sz;
      sx = vxnew[i];
      sy = vynew[i];
      sz = vznew[i];
      vxnew[i] = vxold[i];
      vynew[i] = vyold[i];
      vznew[i] = vzold[i];
      vxold[i] = sx;
      vyold[i] = sy;
      vzold[i] = sz;
   }
}


void propagate_pos_lp_acc(time_prec dt, pos_prec* x_lp, pos_prec* y_lp,
                          pos_prec* z_lp, const vel_prec* vx_lp,
                          const vel_prec* vy_lp, const vel_prec* vz_lp,
                          const pos_prec* xold_lp, const pos_prec* yold_lp,
                          const pos_prec* zold_lp, double scale)
{
   const pos_prec s = scale;
   #pragma acc parallel loop independent async\
               deviceptr(x_lp,y_lp,z_lp,vx_lp,vy_lp,vz_lp,\
               xold_lp,yold_lp,zold_lp)
   for (int i = 0; i < n; ++i) {
      x_lp[i] = xold_lp[i] + vx_lp[i] * dt + s * (x_lp[i] + xold_lp[i]);
      y_lp[i] = yold_lp[i] + vy_lp[i] * dt + s * (y_lp[i] + yold_lp[i]);
      z_lp[i] = zold_lp[i] + vz_lp[i] * dt + s * (z_lp[i] + zold_lp[i]);
   }
}


void propagate_pos_lp2_acc(time_prec dt, const pos_prec* x_lp,
                           const pos_prec* y_lp, const pos_prec* z_lp,
                           pos_prec* xold_lp, pos_prec* yold_lp,
                           pos_prec* zold_lp, double scale)
{
   const pos_prec s = scale;
   #pragma acc parallel loop independent async\
               deviceptr(x_lp,y_lp,z_lp,xold_lp,yold_lp,zold_lp)
   for (int i = 0; i < n; ++i) {
      xold_lp[i] = xold_lp[i] + s * (x_lp[i] + xold_lp[i]);
      yold_lp[i] = yold_lp[i] + s * (y_lp[i] + yold_lp[i]);
      zold_lp[i] = zold_lp[i] + s * (z_lp[i] + zold_lp[i]);
   }
}


void propagate_pos_lf_acc(time_prec dt, pos_prec* qx, pos_prec* qy,
                          pos_prec* qz, const pos_prec* qxold,
                          const pos_prec* qyold, const pos_prec* qzold,
                          const vel_prec* vlx, const vel_prec* vly,
                          const vel_prec* vlz)
{
   #pragma acc parallel loop independent async\
               deviceptr(qx,qy,qz,qxold,qyold,qzold,vlx,vly,vlz)
   for (int i = 0; i < n; ++i) {
      qx[i] = qxold[i] + dt * vlx[i];
      qy[i] = qyold[i] + dt * vly[i];
      qz[i] = qzold[i] + dt * vlz[i];
   }
}


void propagate_velocity_lp_acc(
   vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp, const vel_prec* vxnew_lp,
   const vel_prec* vynew_lp, const vel_prec* vznew_lp, const vel_prec* vxold_lp,
   const vel_prec* vyold_lp, const vel_prec* vzold_lp, const double scale,
   energy_prec& eksum_new, energy_prec& eksum_old)
{
   const energy_prec ekcal_inv = 1.0 / units::ekcal;
   energy_prec ekn = 0, eko = 0;
   #pragma acc parallel loop independent async\
               copy(ekn,eko) reduction(+:ekn,eko)\
               deviceptr(mass,vx_lp,vy_lp,vz_lp,vxnew_lp,vynew_lp,vznew_lp,\
               vxold_lp,vyold_lp,vzold_lp)
   for (int i = 0; i < n; ++i) {
      vel_prec vxf = 0.5f * (vx_lp[i] + vxold_lp[i]);
      vel_prec vyf = 0.5f * (vy_lp[i] + vyold_lp[i]);
      vel_prec vzf = 0.5f * (vz_lp[i] + vzold_lp[i]);
      vx_lp[i] = vxnew_lp[i] + scale * vxf;
      vy_lp[i] = vynew_lp[i] + scale * vyf;
      vz_lp[i] = vznew_lp[i] + scale * vzf;
      energy_prec term = 0.5f * mass[i] * ekcal_inv;
      eko += term * (vxf * vxf + vyf * vyf + vzf * vzf);
      ekn += term *
         (vx_lp[i] * vx_lp[i] + vy_lp[i] * vy_lp[i] + vz_lp[i] * vz_lp[i]);
   }
   #pragma acc wait
   eksum_old = eko;
   eksum_new = ekn;
}


void propagate_velocity_lp2_acc(time_prec dt, vel_prec* vx_lp, vel_prec* vy_lp,
                                vel_prec* vz_lp, const pos_prec* x_lp,
                                const pos_prec* y_lp, const pos_prec* z_lp,
                                const pos_prec* xold_lp,
                                const pos_prec* yold_lp,
                                const pos_prec* zold_lp)
{
   time_prec invdt = 1 / dt;
   #pragma acc parallel loop independent async\
               deviceptr(vx_lp,vy_lp,vz_lp,x_lp,y_lp,z_lp,\
               xold_lp,yold_lp,zold_lp)
   for (int i = 0; i < n; ++i) {
      vx_lp[i] = (x_lp[i] - xold_lp[i]) * invdt;
      vy_lp[i] = (y_lp[i] - yold_lp[i]) * invdt;
      vz_lp[i] = (z_lp[i] - zold_lp[i]) * invdt;
   }
}


void propagate_velocity_lp3_acc(
   vel_prec* vx_lp, vel_prec* vy_lp, vel_prec* vz_lp, const vel_prec* vxnew_lp,
   const vel_prec* vynew_lp, const vel_prec* vznew_lp, const vel_prec* vxold_lp,
   const vel_prec* vyold_lp, const vel_prec* vzold_lp, energy_prec& eksum_new)
{
   const energy_prec ekcal_inv = 1.0 / units::ekcal;
   energy_prec ekn = 0;
   #pragma acc parallel loop independent async\
               copy(ekn) reduction(+:ekn)\
               deviceptr(mass,vx_lp,vy_lp,vz_lp,vxnew_lp,vynew_lp,vznew_lp,\
               vxold_lp,vyold_lp,vzold_lp)
   for (int i = 0; i < n; ++i) {
      vx_lp[i] = 0.5f * (vxnew_lp[i] + vxold_lp[i]);
      vy_lp[i] = 0.5f * (vynew_lp[i] + vyold_lp[i]);
      vz_lp[i] = 0.5f * (vznew_lp[i] + vzold_lp[i]);
      energy_prec term = 0.5f * mass[i] * ekcal_inv;
      ekn += term *
         (vx_lp[i] * vx_lp[i] + vy_lp[i] * vy_lp[i] + vz_lp[i] * vz_lp[i]);
   }
   #pragma acc wait
   eksum_new = ekn;
}
}
