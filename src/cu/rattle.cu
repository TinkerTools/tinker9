#include "add.h"
#include "ff/rattle.h"
#include "launch.h"
#include "md.h"
#include <tinker/detail/units.hh>

namespace tinker {
template <class HTYPE>
__global__
void constrain_methyl_cu1(double eps, int nratch2, const int (*restrict iratch2)[3],
   const pos_prec (*restrict kratch2)[2], int nratch3, const int (*restrict iratch3)[4],
   const pos_prec (*restrict kratch3)[3],

   time_prec dt, pos_prec* restrict xnew, pos_prec* restrict ynew, pos_prec* restrict znew,
   const pos_prec* restrict xold, const pos_prec* restrict yold, const pos_prec* restrict zold,
   const double* restrict massinv, vel_prec* restrict vx, vel_prec* restrict vy,
   vel_prec* restrict vz)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int stride = blockDim.x * gridDim.x;
   constexpr int maxiter = 500;

   int n23 = nratch2 + nratch3;
   for (int im0 = ithread; im0 < n23; im0 += stride) {
      bool methyl = im0 >= nratch2;
      int ia, ib, ic, id;
      double lab, lac, lad;
      double rma, rmb, rmc, rmd;
      if (methyl) {
         int im = im0 - nratch2;
         ia = iratch3[im][0];
         ib = iratch3[im][1];
         ic = iratch3[im][2];
         id = iratch3[im][3];
         lab = kratch3[im][0];
         lac = kratch3[im][1];
         lad = kratch3[im][2];
         rmd = massinv[id];
      } else {
         int im = im0;
         ia = iratch2[im][0];
         ib = iratch2[im][1];
         ic = iratch2[im][2];
         lab = kratch2[im][0];
         lac = kratch2[im][1];
      }
      rma = massinv[ia];
      rmb = massinv[ib];
      rmc = massinv[ic];

      // vectors AB0, AB1;
      double xb0, yb0, zb0, xb1, yb1, zb1;
      xb0 = xold[ib] - xold[ia];
      yb0 = yold[ib] - yold[ia];
      zb0 = zold[ib] - zold[ia];
      xb1 = xnew[ib] - xnew[ia];
      yb1 = ynew[ib] - ynew[ia];
      zb1 = znew[ib] - znew[ia];

      // vectors AC0, AC1;
      double xc0, yc0, zc0, xc1, yc1, zc1;
      xc0 = xold[ic] - xold[ia];
      yc0 = yold[ic] - yold[ia];
      zc0 = zold[ic] - zold[ia];
      xc1 = xnew[ic] - xnew[ia];
      yc1 = ynew[ic] - ynew[ia];
      zc1 = znew[ic] - znew[ia];

      // vectors AD0, AD1
      double xd0, yd0, zd0, xd1, yd1, zd1;
      if (methyl) {
         xd0 = xold[id] - xold[ia];
         yd0 = yold[id] - yold[ia];
         zd0 = zold[id] - zold[ia];
         xd1 = xnew[id] - xnew[ia];
         yd1 = ynew[id] - ynew[ia];
         zd1 = znew[id] - znew[ia];
      }

      double dxa = 0, dya = 0, dza = 0;
      double dxb = 0, dyb = 0, dzb = 0;
      double dxc = 0, dyc = 0, dzc = 0;
      double dxd = 0, dyd = 0, dzd = 0;

      int iter = 0;
      bool done = false;
      while (not done and iter < maxiter) {
         ++iter;
         done = true;
         double x1, y1, z1, dot, delta, term;

         // AB
         x1 = xb1 + dxb - dxa;
         y1 = yb1 + dyb - dya;
         z1 = zb1 + dzb - dza;
         delta = x1 * x1 + y1 * y1 + z1 * z1 - lab * lab;
         if (fabs(delta) > eps) {
            dot = xb0 * x1 + yb0 * y1 + zb0 * z1;
            term = 0.5 * delta / ((rma + rmb) * dot);
            dxa += term * xb0 * rma;
            dya += term * yb0 * rma;
            dza += term * zb0 * rma;
            dxb -= term * xb0 * rmb;
            dyb -= term * yb0 * rmb;
            dzb -= term * zb0 * rmb;
            done = false;
         }

         // AC
         x1 = xc1 + dxc - dxa;
         y1 = yc1 + dyc - dya;
         z1 = zc1 + dzc - dza;
         delta = x1 * x1 + y1 * y1 + z1 * z1 - lac * lac;
         if (fabs(delta) > eps) {
            dot = xc0 * x1 + yc0 * y1 + zc0 * z1;
            term = 0.5 * delta / ((rma + rmc) * dot);
            dxa += term * xc0 * rma;
            dya += term * yc0 * rma;
            dza += term * zc0 * rma;
            dxc -= term * xc0 * rmc;
            dyc -= term * yc0 * rmc;
            dzc -= term * zc0 * rmc;
            done = false;
         }

         // AD
         if (methyl) {
            x1 = xd1 + dxd - dxa;
            y1 = yd1 + dyd - dya;
            z1 = zd1 + dzd - dza;
            delta = x1 * x1 + y1 * y1 + z1 * z1 - lad * lad;
            if (fabs(delta) > eps) {
               dot = xd0 * x1 + yd0 * y1 + zd0 * z1;
               term = 0.5 * delta / ((rma + rmd) * dot);
               dxa += term * xd0 * rma;
               dya += term * yd0 * rma;
               dza += term * zd0 * rma;
               dxd -= term * xd0 * rmd;
               dyd -= term * yd0 * rmd;
               dzd -= term * zd0 * rmd;
               done = false;
            }
         }
      }

      xnew[ia] += dxa;
      ynew[ia] += dya;
      znew[ia] += dza;
      xnew[ib] += dxb;
      ynew[ib] += dyb;
      znew[ib] += dzb;
      xnew[ic] += dxc;
      ynew[ic] += dyc;
      znew[ic] += dzc;
      if (methyl) {
         xnew[id] += dxd;
         ynew[id] += dyd;
         znew[id] += dzd;
      }
      if CONSTEXPR (not eq<HTYPE, SHAKE>()) {
         double invdt = 1 / dt;
         vx[ia] += dxa * invdt;
         vy[ia] += dya * invdt;
         vz[ia] += dza * invdt;
         vx[ib] += dxb * invdt;
         vy[ib] += dyb * invdt;
         vz[ib] += dzb * invdt;
         vx[ic] += dxc * invdt;
         vy[ic] += dyc * invdt;
         vz[ic] += dzc * invdt;
         if (methyl) {
            vx[id] += dxd * invdt;
            vy[id] += dyd * invdt;
            vz[id] += dzd * invdt;
         }
      }
   }
}

void rattle_methyl_cu(
   time_prec dt, const pos_prec* xold, const pos_prec* yold, const pos_prec* zold)
{
   int n23 = nratch2 + nratch3;
   if (n23 <= 0)
      return;

   auto ker = constrain_methyl_cu1<RATTLE>;
   launch_k2s(g::s0, 64, n23, ker,

      rateps, nratch2, iratch2, kratch2, nratch3, iratch3, kratch3,

      dt, xpos, ypos, zpos, xold, yold, zold, massinv, vx, vy, vz);
}

void shake_methyl_cu(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
   const pos_prec* xold, const pos_prec* yold, const pos_prec* zold)
{
   int n23 = nratch2 + nratch3;
   if (n23 <= 0)
      return;

   auto ker = constrain_methyl_cu1<SHAKE>;
   launch_k2s(g::s0, 64, n23, ker,

      rateps, nratch2, iratch2, kratch2, nratch3, iratch3, kratch3,

      dt, xnew, ynew, znew, xold, yold, zold, massinv, nullptr, nullptr, nullptr);
}

template <bool DO_V>
__global__
void constrain2_methyl_cu1(int nratch2, const int (*restrict iratch2)[3], int nratch3,
   const int (*restrict iratch3)[4],

   time_prec dt, vel_prec* restrict vx, vel_prec* restrict vy, vel_prec* restrict vz,
   virial_buffer restrict vir_buf,

   const pos_prec* restrict xpos, const pos_prec* restrict ypos, const pos_prec* restrict zpos,
   const double* restrict massinv)
{
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int stride = blockDim.x * gridDim.x;

   const double vterm = 2 / (dt * units::ekcal);
   double vxx, vyx, vzx, vyy, vzy, vzz;
   if CONSTEXPR (DO_V) {
      vxx = 0, vyx = 0, vzx = 0, vyy = 0, vzy = 0, vzz = 0;
   }

   int n23 = nratch2 + nratch3;
   for (int im0 = ithread; im0 < n23; im0 += stride) {
      bool methyl = im0 >= nratch2;
      int ia, ib, ic, id;
      double rma, rmb, rmc, rmd;
      if (methyl) {
         int im = im0 - nratch2;
         ia = iratch3[im][0];
         ib = iratch3[im][1];
         ic = iratch3[im][2];
         id = iratch3[im][3];
         rmd = massinv[id];
      } else {
         int im = im0;
         ia = iratch2[im][0];
         ib = iratch2[im][1];
         ic = iratch2[im][2];
      }
      rma = massinv[ia];
      rmb = massinv[ib];
      rmc = massinv[ic];

      // matrix form
      // (mab AB3.AB3   rma AB3.AC3  rma AB3.AD3) (lb) = (AB3.vAB0)
      // (rma AC3.AB3   mac AC3.AC3  rma AC3.AD3) (lc) = (AC3.vAC0)
      // (rma AD3.AB3   rma AD3.AC3  mad AD3.AD3) (ld) = (AD3.vAD0)

      // vectors AB3, vAB0, AB3 dot vAB0
      double xb3, yb3, zb3, vxb, vyb, vzb, dotb, rb2;
      xb3 = xpos[ib] - xpos[ia];
      yb3 = ypos[ib] - ypos[ia];
      zb3 = zpos[ib] - zpos[ia];
      vxb = vx[ib] - vx[ia];
      vyb = vy[ib] - vy[ia];
      vzb = vz[ib] - vz[ia];
      dotb = xb3 * vxb + yb3 * vyb + zb3 * vzb;
      rb2 = xb3 * xb3 + yb3 * yb3 + zb3 * zb3;
      // vectors AC3, vAC0, AC3 dot vAC0
      double xc3, yc3, zc3, vxc, vyc, vzc, dotc, rc2;
      xc3 = xpos[ic] - xpos[ia];
      yc3 = ypos[ic] - ypos[ia];
      zc3 = zpos[ic] - zpos[ia];
      vxc = vx[ic] - vx[ia];
      vyc = vy[ic] - vy[ia];
      vzc = vz[ic] - vz[ia];
      dotc = xc3 * vxc + yc3 * vyc + zc3 * vzc;
      rc2 = xc3 * xc3 + yc3 * yc3 + zc3 * zc3;
      // AB3 dot AC3
      double dotbc = xb3 * xc3 + yb3 * yc3 + zb3 * zc3;
      // vectors AD3, vAD0, AD3 dot vAD0
      // AB3 dot AD3, AC3 dot AD3
      double xd3 = 0, yd3 = 0, zd3 = 0, vxd, vyd, vzd, dotd, rd2;
      double dotbd, dotcd;
      if (methyl) {
         xd3 = xpos[id] - xpos[ia];
         yd3 = ypos[id] - ypos[ia];
         zd3 = zpos[id] - zpos[ia];
         vxd = vx[id] - vx[ia];
         vyd = vy[id] - vy[ia];
         vzd = vz[id] - vz[ia];
         dotd = xd3 * vxd + yd3 * vyd + zd3 * vzd;
         rd2 = xd3 * xd3 + yd3 * yd3 + zd3 * zd3;
         dotbd = xb3 * xd3 + yb3 * yd3 + zb3 * zd3;
         dotcd = xc3 * xd3 + yc3 * yd3 + zc3 * zd3;
      }

      double lb, lc, ld;
      double m11, m12, m22; // m21 = m12
      double m13, m23, m33; // m31 = m13, m32 = m23
      double det;
      m11 = (rma + rmb) * rb2;
      m12 = rma * dotbc;
      m22 = (rma + rmc) * rc2;
      if (not methyl) {
         det = m11 * m22 - m12 * m12;
         det = 1 / det;
         lb = (m22 * dotb - m12 * dotc) * det;
         lc = (m11 * dotc - m12 * dotb) * det;
         ld = 0;
      } else {
         m13 = rma * dotbd;
         m23 = rma * dotcd;
         m33 = (rma + rmd) * rd2;
         det = (m11 * m22 - m12 * m12) * m33 + (m12 * m13 - m11 * m23) * m23 +
            (m12 * m23 - m22 * m13) * m13;
         det = 1 / det;
         double i11 = m22 * m33 - m23 * m23;
         double i22 = m11 * m33 - m13 * m13;
         double i33 = m11 * m22 - m12 * m12;
         double i12 = m13 * m23 - m12 * m33;
         double i13 = m12 * m23 - m13 * m22;
         double i23 = m12 * m13 - m11 * m23;
         lb = (i11 * dotb + i12 * dotc + i13 * dotd) * det;
         lc = (i12 * dotb + i22 * dotc + i23 * dotd) * det;
         ld = (i13 * dotb + i23 * dotc + i33 * dotd) * det;
      }

      lb = -lb;
      lc = -lc;
      ld = -ld;
      double xtermb, ytermb, ztermb, xtermc, ytermc, ztermc;
      xtermb = xb3 * lb;
      ytermb = yb3 * lb;
      ztermb = zb3 * lb;
      xtermc = xc3 * lc;
      ytermc = yc3 * lc;
      ztermc = zc3 * lc;
      double xtermd = 0, ytermd = 0, ztermd = 0;
      if (methyl) {
         xtermd = xd3 * ld;
         ytermd = yd3 * ld;
         ztermd = zd3 * ld;
      }
      vx[ia] -= (xtermb + xtermc + xtermd) * rma;
      vy[ia] -= (ytermb + ytermc + ytermd) * rma;
      vz[ia] -= (ztermb + ztermc + ztermd) * rma;
      vx[ib] += xtermb * rmb;
      vy[ib] += ytermb * rmb;
      vz[ib] += ztermb * rmb;
      vx[ic] += xtermc * rmc;
      vy[ic] += ytermc * rmc;
      vz[ic] += ztermc * rmc;
      if (methyl) {
         vx[id] += xtermd * rmd;
         vy[id] += ytermd * rmd;
         vz[id] += ztermd * rmd;
      }
      if CONSTEXPR (DO_V) {
         xtermb *= vterm;
         ytermb *= vterm;
         ztermb *= vterm;
         xtermc *= vterm;
         ytermc *= vterm;
         ztermc *= vterm;
         xtermd *= vterm;
         ytermd *= vterm;
         ztermd *= vterm;
         vxx -= (xb3 * xtermb + xc3 * xtermc + xd3 * xtermd);
         vyx -= (yb3 * xtermb + yc3 * xtermc + yd3 * xtermd);
         vzx -= (zb3 * xtermb + zc3 * xtermc + zd3 * xtermd);
         vyy -= (yb3 * ytermb + yc3 * ytermc + yd3 * ytermd);
         vzz -= (zb3 * ztermb + zc3 * ztermc + zd3 * ztermd);
         vzy -= (zb3 * ytermb + zc3 * ytermc + zd3 * ytermd);
      }
   }

   if CONSTEXPR (DO_V) {
      atomic_add(
         (real)vxx, (real)vyx, (real)vzx, (real)vyy, (real)vzy, (real)vzz, vir_buf, ithread);
   }
}

void rattle2_methyl_cu(time_prec dt, bool do_v)
{
   int n23 = nratch2 + nratch3;
   if (n23 <= 0)
      return;

   if (do_v) {
      auto ker = constrain2_methyl_cu1<true>;
      launch_k2b(g::s0, 64, n23, ker,

         nratch2, iratch2, nratch3, iratch3,

         dt, vx, vy, vz, vir_buf,

         xpos, ypos, zpos, massinv);
   } else {
      auto ker = constrain2_methyl_cu1<false>;
      launch_k2b(g::s0, 64, n23, ker,

         nratch2, iratch2, nratch3, iratch3,

         dt, vx, vy, vz, nullptr,

         xpos, ypos, zpos, massinv);
   }
}
}
