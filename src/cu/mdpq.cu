#include "md/pq.h"
#include "seq/add.h"
#include "seq/launch.h"
#include <tinker/detail/units.hh>

namespace tinker {
__global__
void mdPos_cu1(int n, time_prec dt,                                     //
   pos_prec* restrict qx, pos_prec* restrict qy, pos_prec* restrict qz, //
   const vel_prec* restrict vlx, const vel_prec* restrict vly, const vel_prec* restrict vlz)
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      qx[i] += dt * vlx[i];
      qy[i] += dt * vly[i];
      qz[i] += dt * vlz[i];
   }
}

void mdPos_cu(time_prec dt, pos_prec* qx, pos_prec* qy, pos_prec* qz, const vel_prec* vlx,
   const vel_prec* vly, const vel_prec* vlz)
{
   launch_k1s(g::s0, n, mdPos_cu1, n, dt, qx, qy, qz, vlx, vly, vlz);
}
}

namespace tinker {
__global__
void mdVel_cu1(int n, time_prec dt, const double* restrict massinv,        //
   vel_prec* restrict vlx, vel_prec* restrict vly, vel_prec* restrict vlz, //
   const grad_prec* restrict grx, const grad_prec* restrict gry, const grad_prec* restrict grz)
{
   const vel_prec ekcal = units::ekcal;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      vel_prec coef = -ekcal * massinv[i] * dt;
      auto gxi = toFloatGrad<vel_prec>(grx[i]);
      auto gyi = toFloatGrad<vel_prec>(gry[i]);
      auto gzi = toFloatGrad<vel_prec>(grz[i]);
      vlx[i] += coef * gxi;
      vly[i] += coef * gyi;
      vlz[i] += coef * gzi;
   }
}

void mdVel_cu(time_prec dt, vel_prec* vlx, vel_prec* vly, vel_prec* vlz, const grad_prec* grx,
   const grad_prec* gry, const grad_prec* grz)
{
   launch_k1s(g::s0, n, mdVel_cu1, //
      n, dt, massinv, vlx, vly, vlz, grx, gry, grz);
}

__global__
void mdVel2_cu1(int n, const double* restrict massinv, vel_prec* restrict vlx,
   vel_prec* restrict vly, vel_prec* restrict vlz, //
   time_prec dt1, const grad_prec* restrict grx1, const grad_prec* restrict gry1,
   const grad_prec* restrict grz1, //
   time_prec dt2, const grad_prec* restrict grx2, const grad_prec* restrict gry2,
   const grad_prec* restrict grz2)
{
   const vel_prec ekcal = units::ekcal;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      vel_prec coef = -ekcal * massinv[i];
      auto gx1 = toFloatGrad<vel_prec>(grx1[i]);
      auto gy1 = toFloatGrad<vel_prec>(gry1[i]);
      auto gz1 = toFloatGrad<vel_prec>(grz1[i]);
      auto gx2 = toFloatGrad<vel_prec>(grx2[i]);
      auto gy2 = toFloatGrad<vel_prec>(gry2[i]);
      auto gz2 = toFloatGrad<vel_prec>(grz2[i]);
      vlx[i] += coef * (gx1 * dt1 + gx2 * dt2);
      vly[i] += coef * (gy1 * dt1 + gy2 * dt2);
      vlz[i] += coef * (gz1 * dt1 + gz2 * dt2);
   }
}

void mdVel2_cu(time_prec dt, const grad_prec* grx, const grad_prec* gry, const grad_prec* grz,
   time_prec dt2, const grad_prec* grx2, const grad_prec* gry2, const grad_prec* grz2)
{
   launch_k1s(g::s0, n, mdVel2_cu1, //
      n, massinv, vx, vy, vz,       //
      dt, grx, gry, grz,            //
      dt2, grx2, gry2, grz2);
}
}
