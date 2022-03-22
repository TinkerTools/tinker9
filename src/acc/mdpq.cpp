#include "add.h"
#include "image.h"
#include "md.h"
#include "molecule.h"
#include <tinker/detail/units.hh>

namespace tinker {
void mdCopyPosToXyz_acc()
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

void mdPos_acc(time_prec dt, pos_prec* qx, pos_prec* qy, pos_prec* qz, const vel_prec* vlx,
   const vel_prec* vly, const vel_prec* vlz)
{
   #pragma acc parallel loop independent async\
               deviceptr(qx,qy,qz,vlx,vly,vlz)
   for (int i = 0; i < n; ++i) {
      qx[i] += dt * vlx[i];
      qy[i] += dt * vly[i];
      qz[i] += dt * vlz[i];
   }
}

void propagate_pos_axbv_acc(pos_prec a, pos_prec b)
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

void mdBounds_acc()
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

void mdVelB_acc(time_prec dt, vel_prec* vlx, vel_prec* vly, vel_prec* vlz, const vel_prec* vlx0,
   const vel_prec* vly0, const vel_prec* vlz0, const grad_prec* grx, const grad_prec* gry,
   const grad_prec* grz)
{
   const vel_prec ekcal = units::ekcal;
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vlx,vly,vlz,vlx0,vly0,vlz0,grx,gry,grz)
   for (int i = 0; i < n; ++i) {
      vel_prec coef = -ekcal * massinv[i] * dt;
#if TINKER_DETERMINISTIC_FORCE
      vlx[i] = vlx0[i] + coef * fixedTo<vel_prec>(grx[i]);
      vly[i] = vly0[i] + coef * fixedTo<vel_prec>(gry[i]);
      vlz[i] = vlz0[i] + coef * fixedTo<vel_prec>(grz[i]);
#else
      vlx[i] = vlx0[i] + coef * grx[i];
      vly[i] = vly0[i] + coef * gry[i];
      vlz[i] = vlz0[i] + coef * grz[i];
#endif
   }
}

void mdVelA_acc(time_prec dt, vel_prec* vlx, vel_prec* vly, vel_prec* vlz, const grad_prec* grx,
   const grad_prec* gry, const grad_prec* grz)
{
   const vel_prec ekcal = units::ekcal;
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vlx,vly,vlz,grx,gry,grz)
   for (int i = 0; i < n; ++i) {
      vel_prec coef = -ekcal * massinv[i] * dt;
#if TINKER_DETERMINISTIC_FORCE
      vlx[i] += coef * fixedTo<vel_prec>(grx[i]);
      vly[i] += coef * fixedTo<vel_prec>(gry[i]);
      vlz[i] += coef * fixedTo<vel_prec>(grz[i]);
#else
      vlx[i] += coef * grx[i];
      vly[i] += coef * gry[i];
      vlz[i] += coef * grz[i];
#endif
   }
}

void mdVel2_acc(time_prec dt, const grad_prec* grx, const grad_prec* gry, const grad_prec* grz,
   time_prec dt2, const grad_prec* grx2, const grad_prec* gry2, const grad_prec* grz2)
{
   const vel_prec ekcal = units::ekcal;
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vx,vy,vz,grx,gry,grz,grx2,gry2,grz2)
   for (int i = 0; i < n; ++i) {
      vel_prec coef = -ekcal * massinv[i];
#if TINKER_DETERMINISTIC_FORCE
      // clang-format off
      vx[i] += coef*(fixedTo<vel_prec>(grx[i])*dt+fixedTo<vel_prec>(grx2[i])*dt2);
      vy[i] += coef*(fixedTo<vel_prec>(gry[i])*dt+fixedTo<vel_prec>(gry2[i])*dt2);
      vz[i] += coef*(fixedTo<vel_prec>(grz[i])*dt+fixedTo<vel_prec>(grz2[i])*dt2);
      // clang-format on
#else
      vx[i] += coef * (grx[i] * dt + grx2[i] * dt2);
      vy[i] += coef * (gry[i] * dt + gry2[i] * dt2);
      vz[i] += coef * (grz[i] * dt + grz2[i] * dt2);
#endif
   }
}
}
