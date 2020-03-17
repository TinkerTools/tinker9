#include "image.h"
#include "mdpq.h"
#include "molecule.h"
#include "syntax/acc/add_def.h"
#include <tinker/detail/units.hh>


TINKER_NAMESPACE_BEGIN
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


void bounds_pos_acc()
{
   auto nmol = molecule.nmol;
   const auto* imol = molecule.imol;
   const auto* kmol = molecule.kmol;


   #pragma acc parallel loop independent async\
               deviceptr(imol,kmol,xpos,ypos,zpos)
   for (int i = 0; i < nmol; ++i) {
      // locate the center of each molecule
      pos_prec xmid = 0, ymid = 0, zmid = 0;
      int start = imol[i][0];
      int stop = imol[i][1];
      #pragma acc loop seq
      for (int j = start; j < stop; ++j) {
         int k = kmol[j];
         xmid += xpos[k];
         ymid += ypos[k];
         zmid += zpos[k];
      }
      int weigh = stop - start;
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


void propagate_velocity_acc(time_prec dt, const fixed* grx, const fixed* gry,
                            const fixed* grz)
{
   const vel_prec ekcal = units::ekcal;
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vx,vy,vz,grx,gry,grz)
   for (int i = 0; i < n; ++i) {
      vel_prec coef = -ekcal * massinv[i] * dt;
      vx[i] += coef * to_flt_acc<vel_prec>(grx[i]);
      vy[i] += coef * to_flt_acc<vel_prec>(gry[i]);
      vz[i] += coef * to_flt_acc<vel_prec>(grz[i]);
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


void propagate_velocity2_acc(time_prec dt, const fixed* grx, const fixed* gry,
                             const fixed* grz, time_prec dt2, const fixed* grx2,
                             const fixed* gry2, const fixed* grz2)
{
   const vel_prec ekcal = units::ekcal;
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vx,vy,vz,grx,gry,grz,grx2,gry2,grz2)
   for (int i = 0; i < n; ++i) {
      vel_prec coef = -ekcal * massinv[i];
      // clang-format off
      vx[i] += coef*(to_flt_acc<vel_prec>(grx[i])*dt+to_flt_acc<vel_prec>(grx2[i])*dt2);
      vy[i] += coef*(to_flt_acc<vel_prec>(gry[i])*dt+to_flt_acc<vel_prec>(gry2[i])*dt2);
      vz[i] += coef*(to_flt_acc<vel_prec>(grz[i])*dt+to_flt_acc<vel_prec>(grz2[i])*dt2);
      // clang-format on
   }
}
TINKER_NAMESPACE_END
