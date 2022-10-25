#include "md/pq.h"
#include "seq/add.h"
#include "seq/launch.h"
#include <tinker/detail/units.hh>

namespace tinker {
#include "mdPos_cu1.cc"
void mdPos_cu(time_prec dt, pos_prec* qx, pos_prec* qy, pos_prec* qz, const vel_prec* vlx, const vel_prec* vly,
   const vel_prec* vlz)
{
   launch_k1s(g::s0, n, mdPos_cu1, //
      n, dt, qx, qy, qz, vlx, vly, vlz);
}

__global__
static void mdPosAxbv_cu1(int n, pos_prec sa, pos_prec sb, //
   pos_prec* restrict xpos, pos_prec* restrict ypos, pos_prec* restrict zpos, const vel_prec* restrict vx,
   const vel_prec* restrict vy, const vel_prec* restrict vz)
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      xpos[i] = sa * xpos[i] + sb * vx[i];
      ypos[i] = sa * ypos[i] + sb * vy[i];
      zpos[i] = sa * zpos[i] + sb * vz[i];
   }
}

void mdPosAxbv_cu(pos_prec a, pos_prec b)
{
   launch_k1s(g::s0, n, mdPosAxbv_cu1, //
      n, a, b, xpos, ypos, zpos, vx, vy, vz);
}

__global__
static void mdPosAxbvAn_cu1(int n,     //
   double3 a0, double3 a1, double3 a2, //
   double3 b0, double3 b1, double3 b2, //
   pos_prec* restrict xpos, pos_prec* restrict ypos, pos_prec* restrict zpos, const vel_prec* restrict vx,
   const vel_prec* restrict vy, const vel_prec* restrict vz)
{
   pos_prec a00 = a0.x, a01 = a0.y, a02 = a0.z;
   pos_prec a10 = a1.x, a11 = a1.y, a12 = a1.z;
   pos_prec a20 = a2.x, a21 = a2.y, a22 = a2.z;
   pos_prec b00 = b0.x, b01 = b0.y, b02 = b0.z;
   pos_prec b10 = b1.x, b11 = b1.y, b12 = b1.z;
   pos_prec b20 = b2.x, b21 = b2.y, b22 = b2.z;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      auto o = xpos[i], p = ypos[i], q = zpos[i];
      auto r = vx[i], s = vy[i], t = vz[i];
      xpos[i] = (a00 * o + a01 * p + a02 * q) + (b00 * r + b01 * s + b02 * t);
      ypos[i] = (a10 * o + a11 * p + a12 * q) + (b10 * r + b11 * s + b12 * t);
      zpos[i] = (a20 * o + a21 * p + a22 * q) + (b20 * r + b21 * s + b22 * t);
   }
}

void mdPosAxbvAn_cu(pos_prec (*a)[3], pos_prec (*b)[3])
{
   double3 a0 = make_double3(a[0][0], a[0][1], a[0][2]);
   double3 a1 = make_double3(a[1][0], a[1][1], a[1][2]);
   double3 a2 = make_double3(a[2][0], a[2][1], a[2][2]);
   double3 b0 = make_double3(b[0][0], b[0][1], b[0][2]);
   double3 b1 = make_double3(b[1][0], b[1][1], b[1][2]);
   double3 b2 = make_double3(b[2][0], b[2][1], b[2][2]);
   launch_k1s(g::s0, n, mdPosAxbvAn_cu1, //
      n, a0, a1, a2, b0, b1, b2, xpos, ypos, zpos, vx, vy, vz);
}
}

namespace tinker {
__global__
static void mdVel_cu1(int n, time_prec dt, const double* restrict massinv, //
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

void mdVel_cu(time_prec dt, vel_prec* vlx, vel_prec* vly, vel_prec* vlz, const grad_prec* grx, const grad_prec* gry,
   const grad_prec* grz)
{
   launch_k1s(g::s0, n, mdVel_cu1, //
      n, dt, massinv, vlx, vly, vlz, grx, gry, grz);
}

__global__
static void mdVel2_cu1(int n, const double* restrict massinv, vel_prec* restrict vlx, vel_prec* restrict vly,
   vel_prec* restrict vlz, //
   time_prec dt1, const grad_prec* restrict grx1, const grad_prec* restrict gry1,
   const grad_prec* restrict grz1, //
   time_prec dt2, const grad_prec* restrict grx2, const grad_prec* restrict gry2, const grad_prec* restrict grz2)
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

void mdVel2_cu(time_prec dt, const grad_prec* grx, const grad_prec* gry, const grad_prec* grz, time_prec dt2,
   const grad_prec* grx2, const grad_prec* gry2, const grad_prec* grz2)
{
   launch_k1s(g::s0, n, mdVel2_cu1, //
      n, massinv, vx, vy, vz,       //
      dt, grx, gry, grz,            //
      dt2, grx2, gry2, grz2);
}

__global__
static void mdVelScale_cu1(int nelem, vel_prec sc, vel_prec* vx0, vel_prec* vy0, vel_prec* vz0)
{
   for (int i = ITHREAD; i < nelem; i += STRIDE) {
      vx0[i] *= sc;
      vy0[i] *= sc;
      vz0[i] *= sc;
   }
}

void mdVelScale_cu(vel_prec sc, int nelem, vel_prec* vx0, vel_prec* vy0, vel_prec* vz0)
{
   launch_k1s(g::s0, nelem, mdVelScale_cu1, nelem, sc, vx0, vy0, vz0);
}

__global__
static void mdVelAvbf_cu1(int n, const double* restrict massinv,        //
   vel_prec a, vel_prec b,                                              //
   vel_prec* restrict vx, vel_prec* restrict vy, vel_prec* restrict vz, //
   const grad_prec* restrict gx1, const grad_prec* restrict gy1, const grad_prec* restrict gz1)
{
   auto ekcal = units::ekcal;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      auto coef = -ekcal * massinv[i];
      auto v0x = vx[i];
      auto v0y = vy[i];
      auto v0z = vz[i];
      auto grx = toFloatGrad<vel_prec>(gx1[i]);
      auto gry = toFloatGrad<vel_prec>(gy1[i]);
      auto grz = toFloatGrad<vel_prec>(gz1[i]);
      vx[i] = a * v0x + coef * b * grx;
      vy[i] = a * v0y + coef * b * gry;
      vz[i] = a * v0z + coef * b * grz;
   }
}

__global__
static void mdVelAvbf_cu2(int n, const double* restrict massinv, int nrespa, //
   vel_prec a, vel_prec b,                                                   //
   vel_prec* restrict vx, vel_prec* restrict vy, vel_prec* restrict vz,      //
   const grad_prec* restrict gx1, const grad_prec* restrict gy1, const grad_prec* restrict gz1,
   const grad_prec* restrict gx2, const grad_prec* restrict gy2, const grad_prec* restrict gz2)
{
   auto ekcal = units::ekcal;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      auto coef = -ekcal * massinv[i];
      auto v0x = vx[i];
      auto v0y = vy[i];
      auto v0z = vz[i];
      auto grx = toFloatGrad<vel_prec>(gx1[i]) / nrespa + toFloatGrad<vel_prec>(gx2[i]);
      auto gry = toFloatGrad<vel_prec>(gy1[i]) / nrespa + toFloatGrad<vel_prec>(gy2[i]);
      auto grz = toFloatGrad<vel_prec>(gz1[i]) / nrespa + toFloatGrad<vel_prec>(gz2[i]);
      vx[i] = a * v0x + coef * b * grx;
      vy[i] = a * v0y + coef * b * gry;
      vz[i] = a * v0z + coef * b * grz;
   }
}

void mdVelAvbf_cu(int nrespa, vel_prec a, vel_prec b,                   //
   vel_prec* restrict vx, vel_prec* restrict vy, vel_prec* restrict vz, //
   const grad_prec* gx1, const grad_prec* gy1, const grad_prec* gz1,    //
   const grad_prec* gx2, const grad_prec* gy2, const grad_prec* gz2)
{
   if (nrespa == 1) {
      launch_k1s(g::s0, n, mdVelAvbf_cu1, //
         n, massinv, a, b,                //
         vx, vy, vz, gx1, gy1, gz1);
   } else {
      launch_k1s(g::s0, n, mdVelAvbf_cu2, //
         n, massinv, nrespa, a, b,        //
         vx, vy, vz, gx1, gy1, gz1, gx2, gy2, gz2);
   }
}

__global__
static void mdVelAvbfAn_cu1(int n, const double* restrict massinv,         //
   double3 a0, double3 a1, double3 a2, double3 b0, double3 b1, double3 b2, //
   vel_prec* restrict vx, vel_prec* restrict vy, vel_prec* restrict vz,    //
   const grad_prec* restrict gx1, const grad_prec* restrict gy1, const grad_prec* restrict gz1)
{
   auto ekcal = units::ekcal;
   pos_prec a00 = a0.x, a01 = a0.y, a02 = a0.z;
   pos_prec a10 = a1.x, a11 = a1.y, a12 = a1.z;
   pos_prec a20 = a2.x, a21 = a2.y, a22 = a2.z;
   pos_prec b00 = b0.x, b01 = b0.y, b02 = b0.z;
   pos_prec b10 = b1.x, b11 = b1.y, b12 = b1.z;
   pos_prec b20 = b2.x, b21 = b2.y, b22 = b2.z;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      auto coef = -ekcal * massinv[i];
      auto v0x = vx[i];
      auto v0y = vy[i];
      auto v0z = vz[i];
      auto grx = toFloatGrad<vel_prec>(gx1[i]);
      auto gry = toFloatGrad<vel_prec>(gy1[i]);
      auto grz = toFloatGrad<vel_prec>(gz1[i]);
      // clang-format off
      vx[i] = a00*v0x + a01*v0y + a02*v0z + coef*(b00*grx+b01*gry+b02*grz);
      vy[i] = a10*v0x + a11*v0y + a12*v0z + coef*(b10*grx+b11*gry+b12*grz);
      vz[i] = a20*v0x + a21*v0y + a22*v0z + coef*(b20*grx+b21*gry+b22*grz);
      // clang-format on
   }
}

__global__
static void mdVelAvbfAn_cu2(int n, const double* restrict massinv, int nrespa, //
   double3 a0, double3 a1, double3 a2, double3 b0, double3 b1, double3 b2,     //
   vel_prec* restrict vx, vel_prec* restrict vy, vel_prec* restrict vz,        //
   const grad_prec* restrict gx1, const grad_prec* restrict gy1, const grad_prec* restrict gz1,
   const grad_prec* restrict gx2, const grad_prec* restrict gy2, const grad_prec* restrict gz2)
{
   auto ekcal = units::ekcal;
   pos_prec a00 = a0.x, a01 = a0.y, a02 = a0.z;
   pos_prec a10 = a1.x, a11 = a1.y, a12 = a1.z;
   pos_prec a20 = a2.x, a21 = a2.y, a22 = a2.z;
   pos_prec b00 = b0.x, b01 = b0.y, b02 = b0.z;
   pos_prec b10 = b1.x, b11 = b1.y, b12 = b1.z;
   pos_prec b20 = b2.x, b21 = b2.y, b22 = b2.z;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      auto coef = -ekcal * massinv[i];
      auto v0x = vx[i];
      auto v0y = vy[i];
      auto v0z = vz[i];
      auto grx = toFloatGrad<vel_prec>(gx1[i]) / nrespa + toFloatGrad<vel_prec>(gx2[i]);
      auto gry = toFloatGrad<vel_prec>(gy1[i]) / nrespa + toFloatGrad<vel_prec>(gy2[i]);
      auto grz = toFloatGrad<vel_prec>(gz1[i]) / nrespa + toFloatGrad<vel_prec>(gz2[i]);
      // clang-format off
      vx[i] = a00*v0x + a01*v0y + a02*v0z + coef*(b00*grx+b01*gry+b02*grz);
      vy[i] = a10*v0x + a11*v0y + a12*v0z + coef*(b10*grx+b11*gry+b12*grz);
      vz[i] = a20*v0x + a21*v0y + a22*v0z + coef*(b20*grx+b21*gry+b22*grz);
      // clang-format on
   }
}

void mdVelAvbfAn_cu(int nrespa, vel_prec a[3][3], vel_prec b[3][3],     //
   vel_prec* restrict vx, vel_prec* restrict vy, vel_prec* restrict vz, //
   const grad_prec* gx1, const grad_prec* gy1, const grad_prec* gz1,    //
   const grad_prec* gx2, const grad_prec* gy2, const grad_prec* gz2)
{
   double3 a0 = make_double3(a[0][0], a[0][1], a[0][2]);
   double3 a1 = make_double3(a[1][0], a[1][1], a[1][2]);
   double3 a2 = make_double3(a[2][0], a[2][1], a[2][2]);
   double3 b0 = make_double3(b[0][0], b[0][1], b[0][2]);
   double3 b1 = make_double3(b[1][0], b[1][1], b[1][2]);
   double3 b2 = make_double3(b[2][0], b[2][1], b[2][2]);

   if (nrespa == 1) {
      launch_k1s(g::s0, n, mdVelAvbfAn_cu1,  //
         n, massinv, a0, a1, a2, b0, b1, b2, //
         vx, vy, vz, gx1, gy1, gz1);
   } else {
      launch_k1s(g::s0, n, mdVelAvbfAn_cu2,          //
         n, massinv, nrespa, a0, a1, a2, b0, b1, b2, //
         vx, vy, vz, gx1, gy1, gz1, gx2, gy2, gz2);
   }
}
}
