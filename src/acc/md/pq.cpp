#include "md/pq.h"
#include "add.h"
#include "tool/error.h"
#include <tinker/detail/units.hh>

namespace tinker {
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

void mdPosAxbv_acc(pos_prec a, pos_prec b)
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

void mdDebugPosNorm_acc(pos_prec poseps, time_prec dt, //
   const vel_prec* vlx, const vel_prec* vly, const vel_prec* vlz)
{
   int which = -1;
   pos_prec tol2 = poseps * poseps;
   #pragma acc parallel loop independent async\
               copy(which) reduction(max:which)\
               deviceptr(vlx,vly,vlz)
   for (int i = 0; i < n; ++i) {
      pos_prec x1, y1, z1, norm2;
      x1 = dt * vlx[i];
      y1 = dt * vly[i];
      z1 = dt * vlz[i];
      norm2 = x1 * x1 + y1 * y1 + z1 * z1;
      bool big = norm2 > tol2;
      int flag = big ? i : -1;
      which = which > flag ? which : flag;
   }
   #pragma acc wait

   if (which >= 0) {
      printError();
      TINKER_THROW(format("MD-DEBUG POS UPDATE  --  Atom %d Tried to Move More"
                          " than %.4lf Angstroms",
         which + 1, poseps));
   }
}
}

namespace tinker {
void mdVel_acc(time_prec dt, vel_prec* vlx, vel_prec* vly, vel_prec* vlz, const grad_prec* grx,
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

#pragma acc routine seq
static inline vel_prec
#if TINKER_DETERMINISTIC_FORCE
cvt_grad(fixed val)
{
   return fixedTo<vel_prec>(val);
}
#else
cvt_grad(grad_prec val)
{
   return val;
}
#endif

void mdVelAvbf_acc(int nrespa, vel_prec a, vel_prec b, vel_prec* vx, vel_prec* vy, vel_prec* vz,
   const grad_prec* gx1, const grad_prec* gy1, const grad_prec* gz1, const grad_prec* gx2,
   const grad_prec* gy2, const grad_prec* gz2)
{
   auto ekcal = units::ekcal;
   if (nrespa == 1)
      goto label_nrespa1;
   else
      goto label_nrespa2;

label_nrespa1:
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vx,vy,vz,gx1,gy1,gz1)
   for (int i = 0; i < n; ++i) {
      auto coef = -ekcal * massinv[i];
      auto v0x = vx[i];
      auto v0y = vy[i];
      auto v0z = vz[i];
      auto grx = cvt_grad(gx1[i]);
      auto gry = cvt_grad(gy1[i]);
      auto grz = cvt_grad(gz1[i]);
      vx[i] = a * v0x + coef * b * grx;
      vy[i] = a * v0y + coef * b * gry;
      vz[i] = a * v0z + coef * b * grz;
   }
   return;

label_nrespa2:
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vx,vy,vz,gx1,gy1,gz1,gx2,gy2,gz2)
   for (int i = 0; i < n; ++i) {
      auto coef = -ekcal * massinv[i];
      auto v0x = vx[i];
      auto v0y = vy[i];
      auto v0z = vz[i];
      auto grx = cvt_grad(gx1[i]) / nrespa + cvt_grad(gx2[i]);
      auto gry = cvt_grad(gy1[i]) / nrespa + cvt_grad(gy2[i]);
      auto grz = cvt_grad(gz1[i]) / nrespa + cvt_grad(gz2[i]);
      vx[i] = a * v0x + coef * b * grx;
      vy[i] = a * v0y + coef * b * gry;
      vz[i] = a * v0z + coef * b * grz;
   }
   return;
}

void mdVelAvbfAn_acc(int nrespa, vel_prec a[3][3], vel_prec b[3][3], vel_prec* vx, vel_prec* vy,
   vel_prec* vz, const grad_prec* gx1, const grad_prec* gy1, const grad_prec* gz1,
   const grad_prec* gx2, const grad_prec* gy2, const grad_prec* gz2)
{
   auto ekcal = units::ekcal;
   auto a00 = a[0][0], a01 = a[0][1], a02 = a[0][2];
   auto a10 = a[1][0], a11 = a[1][1], a12 = a[1][2];
   auto a20 = a[2][0], a21 = a[2][1], a22 = a[2][2];
   auto b00 = b[0][0], b01 = b[0][1], b02 = b[0][2];
   auto b10 = b[1][0], b11 = b[1][1], b12 = b[1][2];
   auto b20 = b[2][0], b21 = b[2][1], b22 = b[2][2];

   if (nrespa == 1)
      goto label_nrespa1;
   else
      goto label_nrespa2;

label_nrespa1:
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vx,vy,vz,gx1,gy1,gz1)
   for (int i = 0; i < n; ++i) {
      auto coef = -ekcal * massinv[i];
      auto v0x = vx[i];
      auto v0y = vy[i];
      auto v0z = vz[i];
      auto grx = cvt_grad(gx1[i]);
      auto gry = cvt_grad(gy1[i]);
      auto grz = cvt_grad(gz1[i]);
      // clang-format off
      vx[i] = a00*v0x + a01*v0y + a02*v0z + coef*(b00*grx+b01*gry+b02*grz);
      vy[i] = a10*v0x + a11*v0y + a12*v0z + coef*(b10*grx+b11*gry+b12*grz);
      vz[i] = a20*v0x + a21*v0y + a22*v0z + coef*(b20*grx+b21*gry+b22*grz);
      // clang-foramt on
   }
   return;

label_nrespa2:
   #pragma acc parallel loop independent async\
               deviceptr(massinv,vx,vy,vz,gx1,gy1,gz1,gx2,gy2,gz2)
   for (int i = 0; i < n; ++i) {
      auto coef = -ekcal * massinv[i];
      auto v0x = vx[i];
      auto v0y = vy[i];
      auto v0z = vz[i];
      auto grx = cvt_grad(gx1[i]) / nrespa + cvt_grad(gx2[i]);
      auto gry = cvt_grad(gy1[i]) / nrespa + cvt_grad(gy2[i]);
      auto grz = cvt_grad(gz1[i]) / nrespa + cvt_grad(gz2[i]);
      // clang-format off
      vx[i] = a00*v0x + a01*v0y + a02*v0z + coef*(b00*grx+b01*gry+b02*grz);
      vy[i] = a10*v0x + a11*v0y + a12*v0z + coef*(b10*grx+b11*gry+b12*grz);
      vz[i] = a20*v0x + a21*v0y + a22*v0z + coef*(b20*grx+b21*gry+b22*grz);
      // clang-foramt on
   }
   return;
}
}
