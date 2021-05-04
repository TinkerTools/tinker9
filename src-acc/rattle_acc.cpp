#include "add.h"
#include "lpiston.h"
#include "mdpq.h"
#include "nose.h"
#include "rattle.h"
#include "tinker_rt.h"
#include "tool/darray.h"
#include "tool/error.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/units.hh>


namespace tinker {
namespace {
const int maxiter = 500;
const double sor = 1.25;
}


template <class HTYPE>
void constrain_acc(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
                   const pos_prec* xold, const pos_prec* yold,
                   const pos_prec* zold)
{
   if (nratmol <= 0)
      return;


   const double eps = rateps;


   #pragma acc parallel loop independent async\
           deviceptr(massinv,xnew,ynew,znew,vx,vy,vz,irat,krat,iratmol,\
                     xold,yold,zold)
   for (int im = 0; im < nratmol; ++im) {
      int mbegin = iratmol[im][0];
      int mend = iratmol[im][1];


      bool done = false;
      #pragma acc loop seq
      for (int niter = 0; niter < maxiter and not done; ++niter) {
         done = true;
         for (int i = mbegin; i < mend; ++i) {
            int ia = irat[i][0];
            int ib = irat[i][1];
            double ae1 = 1.0, be1 = 1.0;
            if CONSTEXPR (eq<HTYPE, LPRAT>()) {
               ae1 = lp_rats1;
               be1 = lp_rats1;
            }
            double xr, yr, zr, delta;
            xr = xnew[ib] - xnew[ia];
            yr = ynew[ib] - ynew[ia];
            zr = znew[ib] - znew[ia];
            delta = krat[i] * krat[i] - (xr * xr + yr * yr + zr * zr);
            if (fabs(delta) > eps) {
               done = false;
               double xo, yo, zo, dot, rma, rmb, term, xterm, yterm, zterm;
               xo = xold[ib] - xold[ia];
               yo = yold[ib] - yold[ia];
               zo = zold[ib] - zold[ia];
               dot = xr * xo + yr * yo + zr * zo;
               rma = massinv[ia];
               rmb = massinv[ib];
               term = 0.5 * sor * delta / ((rma + rmb) * dot);
               xterm = xo * term;
               yterm = yo * term;
               zterm = zo * term;
               xnew[ia] -= xterm * rma / ae1;
               ynew[ia] -= yterm * rma / ae1;
               znew[ia] -= zterm * rma / ae1;
               xnew[ib] += xterm * rmb / be1;
               ynew[ib] += yterm * rmb / be1;
               znew[ib] += zterm * rmb / be1;
               if CONSTEXPR (not eq<HTYPE, SHAKE>()) {
                  rma /= dt;
                  rmb /= dt;
                  vx[ia] -= xterm * rma;
                  vy[ia] -= yterm * rma;
                  vz[ia] -= zterm * rma;
                  vx[ib] += xterm * rmb;
                  vy[ib] += yterm * rmb;
                  vz[ib] += zterm * rmb;
               }
            } // end if (delta > eps)
         }
      } // end (maxiter or done)
   }
}


template <class HTYPE>
void settle_acc1(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
                 const pos_prec* xold, const pos_prec* yold,
                 const pos_prec* zold)
{
   if (nratwt <= 0)
      return;


   #pragma acc parallel loop independent async\
           deviceptr(xold,yold,zold,xnew,ynew,znew,vx,vy,vz,mass,\
                     iratwt,kratwt)
   for (int iw = 0; iw < nratwt; ++iw) {
      // atoms a, b, c; lengths ab, ac, bc
      int ia, ib, ic;
      double lab, lac, lbc;
      double m0, m1, m2, invm;
      ia = iratwt[iw][0];
      ib = iratwt[iw][1];
      ic = iratwt[iw][2];
      lab = kratwt[iw][0];
      lac = kratwt[iw][1];
      lbc = kratwt[iw][2];
      m0 = mass[ia];
      m1 = mass[ib];
      m2 = mass[ic];
      double ae1 = 1.0, be1 = 1.0, ce1 = 1.0;
      if CONSTEXPR (eq<HTYPE, LPRAT>()) {
         ae1 = lp_rats1;
         be1 = lp_rats1;
         ce1 = lp_rats1;
         m0 *= ae1;
         m1 *= be1;
         m2 *= ce1;
      }
      invm = 1 / (m0 + m1 + m2);


      // ABC0 is the triangle at t0.
      // ABC1 is the triangle at t0+dt in the absence of constraints.


      // global frame vectors AB0 and AC0
      double xb0, yb0, zb0, xc0, yc0, zc0;
      xb0 = xold[ib] - xold[ia];
      yb0 = yold[ib] - yold[ia];
      zb0 = zold[ib] - zold[ia];
      xc0 = xold[ic] - xold[ia];
      yc0 = yold[ic] - yold[ia];
      zc0 = zold[ic] - zold[ia];


      // global frame centroid D1
      double xcom, ycom, zcom;
      xcom = (m0 * xnew[ia] + m1 * xnew[ib] + m2 * xnew[ic]) * invm;
      ycom = (m0 * ynew[ia] + m1 * ynew[ib] + m2 * ynew[ic]) * invm;
      zcom = (m0 * znew[ia] + m1 * znew[ib] + m2 * znew[ic]) * invm;


      // global frame vectors DA1, DB1, DC1
      double xa1, ya1, za1, xb1, yb1, zb1, xc1, yc1, zc1;
      xa1 = xnew[ia] - xcom;
      ya1 = ynew[ia] - ycom;
      za1 = znew[ia] - zcom;
      xb1 = xnew[ib] - xcom;
      yb1 = ynew[ib] - ycom;
      zb1 = znew[ib] - zcom;
      xc1 = xnew[ic] - xcom;
      yc1 = ynew[ic] - ycom;
      zc1 = znew[ic] - zcom;


      // local frame unit vectors x', y', z' and rotation matrix
      double xakszd, yakszd, zakszd;
      double xaksxd, yaksxd, zaksxd;
      double xaksyd, yaksyd, zaksyd;
      double axlng, aylng, azlng;
      double trns11, trns21, trns31;
      double trns12, trns22, trns32;
      double trns13, trns23, trns33;
      // z' = AB0 cross AC0
      xakszd = yb0 * zc0 - zb0 * yc0;
      yakszd = zb0 * xc0 - xb0 * zc0;
      zakszd = xb0 * yc0 - yb0 * xc0;
      // x' = DA1 cross z'
      xaksxd = ya1 * zakszd - za1 * yakszd;
      yaksxd = za1 * xakszd - xa1 * zakszd;
      zaksxd = xa1 * yakszd - ya1 * xakszd;
      // y' = z' cross x'
      xaksyd = yakszd * zaksxd - zakszd * yaksxd;
      yaksyd = zakszd * xaksxd - xakszd * zaksxd;
      zaksyd = xakszd * yaksxd - yakszd * xaksxd;
      // normalize x', y', z' to get rotation matrix
      axlng = 1 / sqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
      aylng = 1 / sqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
      azlng = 1 / sqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);
      trns11 = xaksxd * axlng;
      trns21 = yaksxd * axlng;
      trns31 = zaksxd * axlng;
      trns12 = xaksyd * aylng;
      trns22 = yaksyd * aylng;
      trns32 = zaksyd * aylng;
      trns13 = xakszd * azlng;
      trns23 = yakszd * azlng;
      trns33 = zakszd * azlng;


      // local frame vectors AB0 and AC0
      double xb0d, yb0d;
      double xc0d, yc0d;
      xb0d = trns11 * xb0 + trns21 * yb0 + trns31 * zb0;
      yb0d = trns12 * xb0 + trns22 * yb0 + trns32 * zb0;
      xc0d = trns11 * xc0 + trns21 * yc0 + trns31 * zc0;
      yc0d = trns12 * xc0 + trns22 * yc0 + trns32 * zc0;


      // local frame vectors DA1, DB1, DC1
      // double xa1d, ya1d;
      double za1d;
      double xb1d, yb1d, zb1d;
      double xc1d, yc1d, zc1d;
      // xa1d = trns11 * xa1 + trns21 * ya1 + trns31 * za1;
      // ya1d = trns12 * xa1 + trns22 * ya1 + trns32 * za1;
      za1d = trns13 * xa1 + trns23 * ya1 + trns33 * za1;
      xb1d = trns11 * xb1 + trns21 * yb1 + trns31 * zb1;
      yb1d = trns12 * xb1 + trns22 * yb1 + trns32 * zb1;
      zb1d = trns13 * xb1 + trns23 * yb1 + trns33 * zb1;
      xc1d = trns11 * xc1 + trns21 * yc1 + trns31 * zc1;
      yc1d = trns12 * xc1 + trns22 * yc1 + trns32 * zc1;
      zc1d = trns13 * xc1 + trns23 * yc1 + trns33 * zc1;


      //====================================================================//
      // analogous to Eq.A1 for arbitrary triangular constraint
      // local frame coordinates of abc0'
      // a0' = (0,yra,0); b0' = (xrb,yrb,0); c0' = (xrc,yrc,0)
      // In Eq.A1, ra, rb, and rc are lengths.
      // Here we use the SIGNED coordinates.
      double yra, xrb, yrb, xrc, yrc;
      {
         // first let A = (0,0), B = (lab,0), C = lac(cosA,sinA)
         // centroid G = (xg,yg)
         double cosA, sinA, xg, yg;
         cosA = (lab * lab + lac * lac - lbc * lbc) / (2 * lab * lac);
         sinA = sqrt(1 - cosA * cosA);
         xg = (m1 * lab + m2 * lac * cosA) * invm;
         yg = m2 * lac * sinA * invm;


         // vectors GA, GB, GC
         double xga, yga, xgb, ygb, xgc, ygc;
         double lga, lgb, lgc;
         double cosAGB, cosAGC, sinAGB, sinAGC;
         xga = -xg;
         yga = -yg;
         xgb = lab - xg;
         ygb = -yg;
         xgc = lac * cosA - xg;
         ygc = lac * sinA - yg;
         lga = sqrt(xga * xga + yga * yga);
         lgb = sqrt(xgb * xgb + ygb * ygb);
         lgc = sqrt(xgc * xgc + ygc * ygc);
         cosAGB = (xga * xgb + yga * ygb) / (lga * lgb);
         cosAGC = (xga * xgc + yga * ygc) / (lga * lgc);
         sinAGB = sqrt(1 - cosAGB * cosAGB);
         sinAGC = sqrt(1 - cosAGC * cosAGC);


         yra = lga;
         xrb = -sinAGB * lgb;
         yrb = cosAGB * lgb;
         xrc = sinAGC * lgc;
         yrc = cosAGC * lgc;
      }


      //====================================================================//
      // Eq.A5
      double sinphi, cosphi, sinpsi, cospsi;
      sinphi = za1d / yra;
      cosphi = sqrt(1 - sinphi * sinphi);
      sinpsi = ((zb1d - zc1d) - (yrb - yrc) * sinphi) / ((xrc - xrb) * cosphi);
      cospsi = sqrt(1 - sinpsi * sinpsi);


      //====================================================================//
      // Eq.A3 local frame vectors da2', db2', dc2'
      double ya2d;
      double xb2d, yb2d;
      double xc2d, yc2d;
      ya2d = yra * cosphi;
      xb2d = xrb * cospsi;
      yb2d = yrb * cosphi + xrb * sinphi * sinpsi;
      xc2d = xrc * cospsi;
      yc2d = yrc * cosphi + xrc * sinphi * sinpsi;
      // adjust numerical error for xb2' and xc2'
      // deltx**2 + ybc**2 + zbc**2 = lbc**2, where ybc**2 + zbc**2 = hh2
      double deltx, hh2;
      hh2 = (yc2d - yb2d) * (yc2d - yb2d) + (zc1d - zb1d) * (zc1d - zb1d);
      // adjusted xbc
      deltx = sqrt(lbc * lbc - hh2);
      // adjustment deltx - (xc - xb)
      deltx = deltx - xc2d + xb2d;
      // m0 xa(=0) + m1 xb + m2 xc = 0
      xb2d -= deltx * m2 / (m1 + m2);
      xc2d += deltx * m1 / (m1 + m2);


      //====================================================================//
      // Eq.A15
      double alpa, beta, gama, al2be2, sinthe, costhe;
      alpa = m1 * (yb2d * yb0d + xb2d * xb0d);
      alpa += m2 * (yc2d * yc0d + xc2d * xc0d);
      beta = m1 * (yb2d * xb0d - xb2d * yb0d);
      beta += m2 * (yc2d * xc0d - xc2d * yc0d);
      gama = m1 * (yb1d * xb0d - xb1d * yb0d);
      gama += m2 * (yc1d * xc0d - xc1d * yc0d);
      al2be2 = alpa * alpa + beta * beta;
      sinthe = (alpa * gama - beta * sqrt(al2be2 - gama * gama)) / al2be2;
      costhe = sqrt(1 - sinthe * sinthe);


      //====================================================================//
      // Eq.A4 local frame vectors da3', db3', dc3'
      double xa3d, ya3d, za3d;
      double xb3d, yb3d, zb3d;
      double xc3d, yc3d, zc3d;
      xa3d = -ya2d * sinthe;
      ya3d = ya2d * costhe;
      za3d = za1d;
      xb3d = xb2d * costhe - yb2d * sinthe;
      yb3d = xb2d * sinthe + yb2d * costhe;
      zb3d = zb1d;
      xc3d = xc2d * costhe - yc2d * sinthe;
      yc3d = xc2d * sinthe + yc2d * costhe;
      zc3d = zc1d;


      // global frame coordinates DA3, DB3, DC3
      double xa3, ya3, za3;
      double xb3, yb3, zb3;
      double xc3, yc3, zc3;
      xa3 = trns11 * xa3d + trns12 * ya3d + trns13 * za3d;
      ya3 = trns21 * xa3d + trns22 * ya3d + trns23 * za3d;
      za3 = trns31 * xa3d + trns32 * ya3d + trns33 * za3d;
      xb3 = trns11 * xb3d + trns12 * yb3d + trns13 * zb3d;
      yb3 = trns21 * xb3d + trns22 * yb3d + trns23 * zb3d;
      zb3 = trns31 * xb3d + trns32 * yb3d + trns33 * zb3d;
      xc3 = trns11 * xc3d + trns12 * yc3d + trns13 * zc3d;
      yc3 = trns21 * xc3d + trns22 * yc3d + trns23 * zc3d;
      zc3 = trns31 * xc3d + trns32 * yc3d + trns33 * zc3d;


      xnew[ia] = xcom + xa3;
      ynew[ia] = ycom + ya3;
      znew[ia] = zcom + za3;
      xnew[ib] = xcom + xb3;
      ynew[ib] = ycom + yb3;
      znew[ib] = zcom + zb3;
      xnew[ic] = xcom + xc3;
      ynew[ic] = ycom + yc3;
      znew[ic] = zcom + zc3;


      //====================================================================//
      // velocity correction due to the position constraints


      // This code is not in the original SETTLE paper, but is necessary for
      // the velocity verlet integrator.
      if CONSTEXPR (not eq<HTYPE, SHAKE>()) {
         double invdt = 1 / dt;
         vx[ia] += (xa3 - xa1) * invdt * ae1;
         vy[ia] += (ya3 - ya1) * invdt * ae1;
         vz[ia] += (za3 - za1) * invdt * ae1;
         vx[ib] += (xb3 - xb1) * invdt * be1;
         vy[ib] += (yb3 - yb1) * invdt * be1;
         vz[ib] += (zb3 - zb1) * invdt * be1;
         vx[ic] += (xc3 - xc1) * invdt * ce1;
         vy[ic] += (yc3 - yc1) * invdt * ce1;
         vz[ic] += (zc3 - zc1) * invdt * ce1;
      }
   }
}


template <class HTYPE>
void constrain_ch_acc1(time_prec dt, pos_prec* xnew, pos_prec* ynew,
                       pos_prec* znew, const pos_prec* xold,
                       const pos_prec* yold, const pos_prec* zold)
{
   if (nratch <= 0)
      return;


   #pragma acc parallel loop independent vector_length(64) async\
           deviceptr(xold,yold,zold,xnew,ynew,znew,vx,vy,vz,massinv,\
                     iratch,kratch)
   for (int im = 0; im < nratch; ++im) {
      int ia, ib;
      double lab;
      double rma, rmb, invm;
      ia = iratch[im][0];
      ib = iratch[im][1];
      lab = kratch[im];
      rma = massinv[ia];
      rmb = massinv[ib];
      double ae1 = 1.0, be1 = 1.0;
      if CONSTEXPR (eq<HTYPE, LPRAT>()) {
         ae1 = lp_rats1;
         be1 = lp_rats1;
         rma /= ae1;
         rmb /= be1;
      }
      invm = rma + rmb;


      // vectors AB0, AB1
      // AB0 dot AB1
      double xb0, yb0, zb0, sqb0, xb1, yb1, zb1, sqb1, dot;
      xb0 = xold[ib] - xold[ia];
      yb0 = yold[ib] - yold[ia];
      zb0 = zold[ib] - zold[ia];
      sqb0 = xb0 * xb0 + yb0 * yb0 + zb0 * zb0;
      xb1 = xnew[ib] - xnew[ia];
      yb1 = ynew[ib] - ynew[ia];
      zb1 = znew[ib] - znew[ia];
      sqb1 = xb1 * xb1 + yb1 * yb1 + zb1 * zb1;
      dot = xb0 * xb1 + yb0 * yb1 + zb0 * zb1;


      // al lambda**2 - 2 * be lambda + ga = 0
      double al, be, ga, sqterm, lam;
      al = sqb0 * invm * invm;
      be = invm * dot;
      ga = sqb1 - lab * lab;
      sqterm = be * be - al * ga;
      lam = (be - sqrt(sqterm)) / al;


      // delta_A and delta_B
      double dxa, dya, dza, dxb, dyb, dzb;
      dxa = lam * xb0 * rma;
      dya = lam * yb0 * rma;
      dza = lam * zb0 * rma;
      dxb = -lam * xb0 * rmb;
      dyb = -lam * yb0 * rmb;
      dzb = -lam * zb0 * rmb;


      xnew[ia] += dxa;
      ynew[ia] += dya;
      znew[ia] += dza;
      xnew[ib] += dxb;
      ynew[ib] += dyb;
      znew[ib] += dzb;


      if CONSTEXPR (not eq<HTYPE, SHAKE>()) {
         double invdt = 1 / dt;
         vx[ia] += dxa * invdt * ae1;
         vy[ia] += dya * invdt * ae1;
         vz[ia] += dza * invdt * ae1;
         vx[ib] += dxb * invdt * be1;
         vy[ib] += dyb * invdt * be1;
         vz[ib] += dzb * invdt * be1;
      }
   }
}


void rattle_acc(time_prec dt, const pos_prec* xold, const pos_prec* yold,
                const pos_prec* zold)
{
   constrain_acc<RATTLE>(dt, xpos, ypos, zpos, xold, yold, zold);
}


void rattle_settle_acc(time_prec dt, const pos_prec* xold, const pos_prec* yold,
                       const pos_prec* zold)
{
   settle_acc1<RATTLE>(dt, xpos, ypos, zpos, xold, yold, zold);
}


void rattle_ch_acc(time_prec dt, const pos_prec* xold, const pos_prec* yold,
                   const pos_prec* zold)
{
   constrain_ch_acc1<RATTLE>(dt, xpos, ypos, zpos, xold, yold, zold);
}


void lprat_acc(time_prec dt, const pos_prec* xold, const pos_prec* yold,
               const pos_prec* zold)
{
   constrain_acc<LPRAT>(dt, xpos, ypos, zpos, xold, yold, zold);
}


void lprat_settle_acc(time_prec dt, const pos_prec* xold, const pos_prec* yold,
                      const pos_prec* zold)
{
   settle_acc1<LPRAT>(dt, xpos, ypos, zpos, xold, yold, zold);
}


void lprat_ch_acc(time_prec dt, const pos_prec* xold, const pos_prec* yold,
                  const pos_prec* zold)
{
   constrain_ch_acc1<LPRAT>(dt, xpos, ypos, zpos, xold, yold, zold);
}


void shake_acc(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
               const pos_prec* xold, const pos_prec* yold, const pos_prec* zold)
{
   constrain_acc<SHAKE>(dt, xnew, ynew, znew, xold, yold, zold);
}


void shake_settle_acc(time_prec dt, pos_prec* xnew, pos_prec* ynew,
                      pos_prec* znew, const pos_prec* xold,
                      const pos_prec* yold, const pos_prec* zold)
{
   settle_acc1<SHAKE>(dt, xnew, ynew, znew, xold, yold, zold);
}


void shake_ch_acc(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
                  const pos_prec* xold, const pos_prec* yold,
                  const pos_prec* zold)
{
   constrain_ch_acc1<SHAKE>(dt, xnew, ynew, znew, xold, yold, zold);
}


//====================================================================//


template <class HTYPE, bool DO_V>
void constrain2_acc(time_prec dt)
{
   if (nratmol <= 0)
      return;


   const double eps = rateps / dt;
   const double vterm = 2 / (dt * units::ekcal);
   size_t bufsize = buffer_size();


   #pragma acc parallel loop independent async\
           deviceptr(massinv,xpos,ypos,zpos,vx,vy,vz,vir_buf,irat,krat,iratmol)
   for (int im = 0; im < nratmol; ++im) {
      int mbegin = iratmol[im][0];
      int mend = iratmol[im][1];


      bool done = false;
      #pragma acc loop seq
      for (int niter = 0; niter < maxiter and not done; ++niter) {
         done = true;
         for (int i = mbegin; i < mend; ++i) {
            int ia = irat[i][0];
            int ib = irat[i][1];
            double xr, yr, zr, xv, yv, zv, dot, rma, rmb, term;
            xr = xpos[ib] - xpos[ia];
            yr = ypos[ib] - ypos[ia];
            zr = zpos[ib] - zpos[ia];
            xv = vx[ib] - vx[ia];
            yv = vy[ib] - vy[ia];
            zv = vz[ib] - vz[ia];
            dot = xr * xv + yr * yv + zr * zv;
            rma = massinv[ia];
            rmb = massinv[ib];
            if CONSTEXPR (eq<HTYPE, LPRAT>()) {
               double ae2 = lp_rats2;
               double be2 = lp_rats2;
               rma /= ae2;
               rmb /= be2;
            }
            term = -dot / ((rma + rmb) * krat[i] * krat[i]);
            if (fabs(term) > eps) {
               done = false;
               term *= sor;
               double xterm, yterm, zterm;
               xterm = xr * term;
               yterm = yr * term;
               zterm = zr * term;
               vx[ia] -= xterm * rma;
               vy[ia] -= yterm * rma;
               vz[ia] -= zterm * rma;
               vx[ib] += xterm * rmb;
               vy[ib] += yterm * rmb;
               vz[ib] += zterm * rmb;
               if CONSTEXPR (DO_V) {
                  size_t offset = im & (bufsize - 1);
                  double vxx = 0, vyx = 0, vzx = 0;
                  double vyy = 0, vzy = 0, vzz = 0;
                  xterm *= vterm;
                  yterm *= vterm;
                  zterm *= vterm;
                  vxx -= xr * xterm;
                  vyx -= yr * xterm;
                  vzx -= zr * xterm;
                  vyy -= yr * yterm;
                  vzy -= zr * yterm;
                  vzz -= zr * zterm;
                  atomic_add((real)vxx, (real)vyx, (real)vzx, (real)vyy,
                             (real)vzy, (real)vzz, vir_buf, offset);
               }
            } // end if (delta > eps)
         }
      } // end (maxiter or done)
   }
}


template <class HTYPE, bool DO_V>
void settle2_acc1(time_prec dt)
{
   if (nratwt <= 0)
      return;


   const double vterm = 2 / (dt * units::ekcal);
   size_t bufsize = buffer_size();


   #pragma acc parallel loop independent async\
           deviceptr(vx,vy,vz,xpos,ypos,zpos,mass,vir_buf,iratwt)
   for (int iw = 0; iw < nratwt; ++iw) {
      int ia, ib, ic;
      double m0, m1, m2;
      ia = iratwt[iw][0];
      ib = iratwt[iw][1];
      ic = iratwt[iw][2];
      m0 = mass[ia];
      m1 = mass[ib];
      m2 = mass[ic];
      if CONSTEXPR (eq<HTYPE, LPRAT>()) {
         double ae2 = lp_rats2;
         double be2 = lp_rats2;
         double ce2 = lp_rats2;
         m0 *= ae2;
         m1 *= be2;
         m2 *= ce2;
      }


      // vectors AB, BC, CA
      double xab, yab, zab;
      double xbc, ybc, zbc;
      double xca, yca, zca;
      xab = xpos[ib] - xpos[ia];
      yab = ypos[ib] - ypos[ia];
      zab = zpos[ib] - zpos[ia];
      xbc = xpos[ic] - xpos[ib];
      ybc = ypos[ic] - ypos[ib];
      zbc = zpos[ic] - zpos[ib];
      xca = xpos[ia] - xpos[ic];
      yca = ypos[ia] - ypos[ic];
      zca = zpos[ia] - zpos[ic];


      // unit vectors eAB, eBC, eCA
      double ablng, bclng, calng;
      double xeab, yeab, zeab;
      double xebc, yebc, zebc;
      double xeca, yeca, zeca;
      ablng = 1 / sqrt(xab * xab + yab * yab + zab * zab);
      bclng = 1 / sqrt(xbc * xbc + ybc * ybc + zbc * zbc);
      calng = 1 / sqrt(xca * xca + yca * yca + zca * zca);
      xeab = xab * ablng;
      yeab = yab * ablng;
      zeab = zab * ablng;
      xebc = xbc * bclng;
      yebc = ybc * bclng;
      zebc = zbc * bclng;
      xeca = xca * calng;
      yeca = yca * calng;
      zeca = zca * calng;


      // velocity vectors vAB, vBC, vCA
      double xvab, yvab, zvab;
      double xvbc, yvbc, zvbc;
      double xvca, yvca, zvca;
      xvab = vx[ib] - vx[ia];
      yvab = vy[ib] - vy[ia];
      zvab = vz[ib] - vz[ia];
      xvbc = vx[ic] - vx[ib];
      yvbc = vy[ic] - vy[ib];
      zvbc = vz[ic] - vz[ib];
      xvca = vx[ia] - vx[ic];
      yvca = vy[ia] - vy[ic];
      zvca = vz[ia] - vz[ic];


      double vabab, vbcbc, vcaca;
      vabab = xvab * xeab + yvab * yeab + zvab * zeab;
      vbcbc = xvbc * xebc + yvbc * yebc + zvbc * zebc;
      vcaca = xvca * xeca + yvca * yeca + zvca * zeca;


      double cosa, cosb, cosc;
      cosa = -xeab * xeca - yeab * yeca - zeab * zeca;
      cosb = -xebc * xeab - yebc * yeab - zebc * zeab;
      cosc = -xeca * xebc - yeca * yebc - zeca * zebc;


      //====================================================================//
      // Eq.B1
      // M (tab,tbc,tca) = m2v, where
      // m2v = (ma mb vab0, mb mc vbc0, mc ma vca0),
      //     (a1,a2,a3)
      // M = (b1,b2,b3)
      //     (c1,c2,c3)
      double a1, a2, a3, b1, b2, b3, c1, c2, c3, denom;
      a1 = m0 + m1;
      a2 = m0 * cosb;
      a3 = m1 * cosa;
      b1 = m2 * cosb;
      b2 = m1 + m2;
      b3 = m1 * cosc;
      c1 = m2 * cosa;
      c2 = m0 * cosc;
      c3 = m2 + m0;
      // det(M)
      denom = a1 * (b2 * c3 - b3 * c2) + a2 * (b3 * c1 - b1 * c3) +
         a3 * (b1 * c2 - b2 * c1);


      // inverse(M)*det(M)
      double av1, av2, av3, bv1, bv2, bv3, cv1, cv2, cv3;
      av1 = b2 * c3 - b3 * c2;
      av2 = c2 * a3 - c3 * a2;
      av3 = a2 * b3 - a3 * b2;
      bv1 = b3 * c1 - b1 * c3;
      bv2 = c3 * a1 - c1 * a3;
      bv3 = a3 * b1 - a1 * b3;
      cv1 = b1 * c2 - b2 * c1;
      cv2 = c1 * a2 - c2 * a1;
      cv3 = a1 * b2 - a2 * b1;


      // t = inverse(M)*m2v
      // clang-format off
      double tabd, tbcd, tcad;
      tabd = av1*m0*m1*vabab + av2*m1*m2*vbcbc + av3*m2*m0*vcaca;
      tbcd = bv1*m0*m1*vabab + bv2*m1*m2*vbcbc + bv3*m2*m0*vcaca;
      tcad = cv1*m0*m1*vabab + cv2*m1*m2*vbcbc + cv3*m2*m0*vcaca;
      // clang-format on


      denom = 1 / denom;
      vx[ia] += (xeab * tabd - xeca * tcad) / m0 * denom;
      vz[ia] += (zeab * tabd - zeca * tcad) / m0 * denom;
      vy[ia] += (yeab * tabd - yeca * tcad) / m0 * denom;
      vx[ib] += (xebc * tbcd - xeab * tabd) / m1 * denom;
      vy[ib] += (yebc * tbcd - yeab * tabd) / m1 * denom;
      vz[ib] += (zebc * tbcd - zeab * tabd) / m1 * denom;
      vx[ic] += (xeca * tcad - xebc * tbcd) / m2 * denom;
      vy[ic] += (yeca * tcad - yebc * tbcd) / m2 * denom;
      vz[ic] += (zeca * tcad - zebc * tbcd) / m2 * denom;


      if CONSTEXPR (DO_V) {
         double xterm, yterm, zterm;
         double vxx = 0, vyx = 0, vzx = 0, vyy = 0, vzy = 0, vzz = 0;
         xterm = xeab * tabd * denom * vterm;
         yterm = yeab * tabd * denom * vterm;
         zterm = zeab * tabd * denom * vterm;
         vxx += xab * xterm;
         vyx += yab * xterm;
         vzx += zab * xterm;
         vyy += yab * yterm;
         vzy += zab * yterm;
         vzz += zab * zterm;
         xterm = xebc * tbcd * denom * vterm;
         yterm = yebc * tbcd * denom * vterm;
         zterm = zebc * tbcd * denom * vterm;
         vxx += xbc * xterm;
         vyx += ybc * xterm;
         vzx += zbc * xterm;
         vyy += ybc * yterm;
         vzy += zbc * yterm;
         vzz += zbc * zterm;
         xterm = xeca * tcad * denom * vterm;
         yterm = yeca * tcad * denom * vterm;
         zterm = zeca * tcad * denom * vterm;
         vxx += xca * xterm;
         vyx += yca * xterm;
         vzx += zca * xterm;
         vyy += yca * yterm;
         vzy += zca * yterm;
         vzz += zca * zterm;


         size_t offset = iw & (bufsize - 1);
         atomic_add((real)vxx, (real)vyx, (real)vzx, (real)vyy, (real)vzy,
                    (real)vzz, vir_buf, offset);
      }
   }
}


template <class HTYPE, bool DO_V>
void constrain2_ch_acc1(time_prec dt)
{
   if (nratch <= 0)
      return;


   const double vterm = 2 / (dt * units::ekcal);
   size_t bufsize = buffer_size();


   #pragma acc parallel loop independent vector_length(64) async\
           deviceptr(massinv,xpos,ypos,zpos,vx,vy,vz,vir_buf,iratch)
   for (int im = 0; im < nratch; ++im) {
      int ia, ib;
      double rma, rmb;
      ia = iratch[im][0];
      ib = iratch[im][1];
      rma = massinv[ia];
      rmb = massinv[ib];
      if CONSTEXPR (eq<HTYPE, LPRAT>()) {
         double ae2 = lp_rats2;
         double be2 = lp_rats2;
         rma /= ae2;
         rmb /= be2;
      }


      // vectors AB3, vAB0, AB3 dot vAB0, -lam
      double xr, yr, zr, xv, yv, zv, dot, term;
      xr = xpos[ib] - xpos[ia];
      yr = ypos[ib] - ypos[ia];
      zr = zpos[ib] - zpos[ia];
      xv = vx[ib] - vx[ia];
      yv = vy[ib] - vy[ia];
      zv = vz[ib] - vz[ia];
      dot = xr * xv + yr * yv + zr * zv;
      double r2 = xr * xr + yr * yr + zr * zr;
      term = -dot / ((rma + rmb) * r2);


      double xterm, yterm, zterm;
      xterm = xr * term;
      yterm = yr * term;
      zterm = zr * term;
      vx[ia] -= xterm * rma;
      vy[ia] -= yterm * rma;
      vz[ia] -= zterm * rma;
      vx[ib] += xterm * rmb;
      vy[ib] += yterm * rmb;
      vz[ib] += zterm * rmb;
      if CONSTEXPR (DO_V) {
         size_t offset = im & (bufsize - 1);
         double vxx = 0, vyx = 0, vzx = 0;
         double vyy = 0, vzy = 0, vzz = 0;
         xterm *= vterm;
         yterm *= vterm;
         zterm *= vterm;
         vxx -= xr * xterm;
         vyx -= yr * xterm;
         vzx -= zr * xterm;
         vyy -= yr * yterm;
         vzy -= zr * yterm;
         vzz -= zr * zterm;
         atomic_add((real)vxx, (real)vyx, (real)vzx, (real)vyy, (real)vzy,
                    (real)vzz, vir_buf, offset);
      }
   }
}


void rattle2_acc(time_prec dt, bool do_v)
{
   if (do_v) {
      constrain2_acc<RATTLE, true>(dt);
   } else {
      constrain2_acc<RATTLE, false>(dt);
   }
}


void rattle2_settle_acc(time_prec dt, bool do_v)
{
   if (do_v)
      settle2_acc1<RATTLE, true>(dt);
   else
      settle2_acc1<RATTLE, false>(dt);
}


void rattle2_ch_acc(time_prec dt, bool do_v)
{
   if (do_v) {
      constrain2_ch_acc1<RATTLE, true>(dt);
   } else {
      constrain2_ch_acc1<RATTLE, false>(dt);
   }
}


void lprat2_acc(time_prec dt)
{
   constrain2_acc<LPRAT, false>(dt);
}


void lprat2_settle_acc(time_prec dt)
{
   settle2_acc1<LPRAT, false>(dt);
}


void lprat2_ch_acc(time_prec dt)
{
   constrain2_ch_acc1<LPRAT, false>(dt);
}
}
