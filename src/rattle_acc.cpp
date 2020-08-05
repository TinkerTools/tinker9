#include "mdpq.h"
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


void rattle_acc(time_prec dt, const pos_prec* xold, const pos_prec* yold,
                const pos_prec* zold)
{
   if (nratmol <= 0)
      return;


   int* moved = rattle_moved;
   int* update = rattle_update;
   int* bigeps = rattle_bigdelta;


   const double eps = rateps;
   int niter = 0;
   bool done = false;


   #pragma acc parallel loop independent async\
           deviceptr(moved,update)
   for (int i = 0; i < n; ++i) {
      moved[i] = true;
      update[i] = false;
   }


   while (not done and niter < maxiter) {
      niter += 1;
      done = true;
      #pragma acc parallel loop independent async\
              deviceptr(massinv,xpos,ypos,zpos,vx,vy,vz,irat,krat,iratmol,\
                        moved,update,bigeps,xold,yold,zold)
      for (int im = 0; im < nratmol; ++im) {
         int mbegin = iratmol[im][0];
         int mend = iratmol[im][1];
         #pragma acc loop seq
         for (int i = mbegin; i < mend; ++i) {
            int bigdelta = false;
            int ia = irat[i][0];
            int ib = irat[i][1];
            if (moved[ia] or moved[ib]) {
               pos_prec xr = xpos[ib] - xpos[ia];
               pos_prec yr = ypos[ib] - ypos[ia];
               pos_prec zr = zpos[ib] - zpos[ia];
               pos_prec delta =
                  krat[i] * krat[i] - (xr * xr + yr * yr + zr * zr);
               if (fabs(delta) > eps) {
                  bigdelta = true;
                  update[ia] = true;
                  update[ib] = true;
                  pos_prec xo = xold[ib] - xold[ia];
                  pos_prec yo = yold[ib] - yold[ia];
                  pos_prec zo = zold[ib] - zold[ia];
                  pos_prec dot = xr * xo + yr * yo + zr * zo;
                  mass_prec rma = massinv[ia];
                  mass_prec rmb = massinv[ib];
                  pos_prec term = 0.5f * sor * delta / ((rma + rmb) * dot);
                  pos_prec xterm = xo * term;
                  pos_prec yterm = yo * term;
                  pos_prec zterm = zo * term;
                  xpos[ia] -= xterm * rma;
                  ypos[ia] -= yterm * rma;
                  zpos[ia] -= zterm * rma;
                  xpos[ib] += xterm * rmb;
                  ypos[ib] += yterm * rmb;
                  zpos[ib] += zterm * rmb;
                  rma /= dt;
                  rmb /= dt;
                  vx[ia] -= xterm * rma;
                  vy[ia] -= yterm * rma;
                  vz[ia] -= zterm * rma;
                  vx[ib] += xterm * rmb;
                  vy[ib] += yterm * rmb;
                  vz[ib] += zterm * rmb;
               } // end if (delta > eps)
            }
            bigeps[i] = bigdelta;
         }
      }


      darray::copy(PROCEED_NEW_Q, n, moved, update);
      darray::zero(PROCEED_NEW_Q, n, update);
      int next_iter = parallel::reduce_logic_or(bigeps, nrat, WAIT_NEW_Q);
      done = not next_iter;
   }


   if (niter == maxiter) {
      darray::copy(PROCEED_NEW_Q, n, xpos, xold);
      darray::copy(PROCEED_NEW_Q, n, ypos, yold);
      darray::copy(PROCEED_NEW_Q, n, zpos, zold);
      t_prterr();
      TINKER_THROW("RATTLE  --  Warning, Distance Constraints not Satisfied");
   } else if (inform::debug) {
      print(stdout,
            " RATTLE   --  Distance Constraints met at %5d Iterations\n",
            niter);
   }
}


template <bool DO_V>
void rattle2_acc1(time_prec dt)
{
   if (nratmol <= 0)
      return;


   int* moved = rattle_moved;
   int* update = rattle_update;
   int* bigeps = rattle_bigdelta;


   const double eps = rateps / dt;
   int niter = 0;
   bool done = false;
   const double vterm = 2 / (dt * units::ekcal);


   #pragma acc parallel loop independent async\
           deviceptr(moved,update)
   for (int i = 0; i < n; ++i) {
      moved[i] = true;
      update[i] = false;
   }


   real vxx = 0, vyx = 0, vzx = 0, vyy = 0, vzy = 0, vzz = 0;
   while (not done and niter < maxiter) {
      niter += 1;
      done = true;
      #pragma acc parallel loop independent async\
              deviceptr(massinv,xpos,ypos,zpos,vx,vy,vz,\
                        irat,krat,iratmol,moved,update,bigeps)
      for (int im = 0; im < nratmol; ++im) {
         int mbegin = iratmol[im][0];
         int mend = iratmol[im][1];
         #pragma acc loop seq
         for (int i = mbegin; i < mend; ++i) {
            int bigdelta = false;
            int ia = irat[i][0];
            int ib = irat[i][1];
            if (moved[ia] or moved[ib]) {
               pos_prec xr = xpos[ib] - xpos[ia];
               pos_prec yr = ypos[ib] - ypos[ia];
               pos_prec zr = zpos[ib] - zpos[ia];
               vel_prec xv = vx[ib] - vx[ia];
               vel_prec yv = vy[ib] - vy[ia];
               vel_prec zv = vz[ib] - vz[ia];
               pos_prec dot = xr * xv + yr * yv + zr * zv;
               mass_prec rma = massinv[ia];
               mass_prec rmb = massinv[ib];
               pos_prec term = -dot / ((rma + rmb) * krat[i] * krat[i]);
               if (fabs(term) > eps) {
                  bigdelta = true;
                  update[ia] = true;
                  update[ib] = true;
                  term *= sor;
                  pos_prec xterm = xr * term;
                  pos_prec yterm = yr * term;
                  pos_prec zterm = zr * term;
                  vx[ia] -= xterm * rma;
                  vy[ia] -= yterm * rma;
                  vz[ia] -= zterm * rma;
                  vx[ib] += xterm * rmb;
                  vy[ib] += yterm * rmb;
                  vz[ib] += zterm * rmb;
                  if CONSTEXPR (DO_V) {
                     xterm *= vterm;
                     yterm *= vterm;
                     zterm *= vterm;
                     vxx -= xr * xterm;
                     vyx -= yr * xterm;
                     vzx -= zr * xterm;
                     vyy -= yr * yterm;
                     vzy -= zr * yterm;
                     vzz -= zr * zterm;
                  }
               } // end if (delta > eps)
            }
            bigeps[i] = bigdelta;
         } // end loop nrat
      }


      darray::copy(PROCEED_NEW_Q, n, moved, update);
      darray::zero(PROCEED_NEW_Q, n, update);
      int next_iter = parallel::reduce_logic_or(bigeps, nrat, WAIT_NEW_Q);
      done = not next_iter;
   }


   if CONSTEXPR (DO_V) {
      vir[0] += vxx;
      vir[1] += vyx;
      vir[2] += vzx;
      vir[3] += vyx;
      vir[4] += vyy;
      vir[5] += vzy;
      vir[6] += vzx;
      vir[7] += vzy;
      vir[8] += vzz;
   }


   if (niter == maxiter) {
      t_prterr();
      TINKER_THROW("RATTLE2  --  Warning, Velocity Constraints not Satisfied");
   } else if (inform::debug) {
      print(stdout,
            " RATTLE2  --  Velocity Constraints met at %5d Iterations\n",
            niter);
   }
}


void rattle2_acc(time_prec dt, bool do_v)
{
   if (do_v) {
      rattle2_acc1<true>(dt);
   } else {
      rattle2_acc1<false>(dt);
   }
}


void rattle_settle_acc(time_prec dt, const pos_prec* xold, const pos_prec* yold,
                       const pos_prec* zold)
{
   if (nratwt <= 0)
      return;


   time_prec invdt = 1 / dt;


   #pragma acc parallel loop independent async\
           deviceptr(xold,yold,zold,xpos,ypos,zpos,vx,vy,vz,mass,\
                     iratwt,kratwt)
   for (int iw = 0; iw < nratwt; ++iw) {
      // atoms a, b, c; lengths ab, ac, bc
      int ia, ib, ic;
      pos_prec lab, lac, lbc;
      mass_prec m0, m1, m2, invm;
      ia = iratwt[iw][0];
      ib = iratwt[iw][1];
      ic = iratwt[iw][2];
      lab = kratwt[iw][0];
      lac = kratwt[iw][1];
      lbc = kratwt[iw][2];
      m0 = mass[ia];
      m1 = mass[ib];
      m2 = mass[ic];
      invm = 1 / (m0 + m1 + m2);


      // ABC0 is the triangle at t0.
      // ABC1 is the triangle at t0+dt in the absence of constraints.
      // D1 is the centroid of ABC_1.


      // global frame vectors AB0 and AC0
      pos_prec xb0, yb0, zb0, xc0, yc0, zc0;
      xb0 = xold[ib] - xold[ia];
      yb0 = yold[ib] - yold[ia];
      zb0 = zold[ib] - zold[ia];
      xc0 = xold[ic] - xold[ia];
      yc0 = yold[ic] - yold[ia];
      zc0 = zold[ic] - zold[ia];


      // global frame centroid D1
      pos_prec xcom, ycom, zcom;
      xcom = (m0 * xpos[ia] + m1 * xpos[ib] + m2 * xpos[ic]) * invm;
      ycom = (m0 * ypos[ia] + m1 * ypos[ib] + m2 * ypos[ic]) * invm;
      zcom = (m0 * zpos[ia] + m1 * zpos[ib] + m2 * zpos[ic]) * invm;


      // global frame vectors DA1, DB1, DC1
      pos_prec xa1, ya1, za1, xb1, yb1, zb1, xc1, yc1, zc1;
      xa1 = xpos[ia] - xcom;
      ya1 = ypos[ia] - ycom;
      za1 = zpos[ia] - zcom;
      xb1 = xpos[ib] - xcom;
      yb1 = ypos[ib] - ycom;
      zb1 = zpos[ib] - zcom;
      xc1 = xpos[ic] - xcom;
      yc1 = ypos[ic] - ycom;
      zc1 = zpos[ic] - zcom;


      // local frame unit vectors x', y', z' and rotation matrix
      pos_prec xakszd, yakszd, zakszd;
      pos_prec xaksxd, yaksxd, zaksxd;
      pos_prec xaksyd, yaksyd, zaksyd;
      pos_prec axlng, aylng, azlng;
      pos_prec trns11, trns21, trns31;
      pos_prec trns12, trns22, trns32;
      pos_prec trns13, trns23, trns33;
      // z' = AB_0 cross AC_0
      xakszd = yb0 * zc0 - zb0 * yc0;
      yakszd = zb0 * xc0 - xb0 * zc0;
      zakszd = xb0 * yc0 - yb0 * xc0;
      // x' = DA_1 cross z'
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
      pos_prec xb0d, yb0d;
      pos_prec xc0d, yc0d;
      xb0d = trns11 * xb0 + trns21 * yb0 + trns31 * zb0;
      yb0d = trns12 * xb0 + trns22 * yb0 + trns32 * zb0;
      xc0d = trns11 * xc0 + trns21 * yc0 + trns31 * zc0;
      yc0d = trns12 * xc0 + trns22 * yc0 + trns32 * zc0;


      // local frame vectors DA1, DB1, DC1
      // pos_prec xa1d, ya1d;
      pos_prec za1d;
      pos_prec xb1d, yb1d, zb1d;
      pos_prec xc1d, yc1d, zc1d;
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
      pos_prec yra, xrb, yrb, xrc, yrc;
      {
         // first let A = (0,0), B = (lab,0), C = lac(cosA,sinA)
         // centroid G = (xg,yg)
         pos_prec cosA, sinA, xg, yg;
         cosA = (lab * lab + lac * lac - lbc * lbc) / (2 * lab * lac);
         sinA = sqrt(1 - cosA * cosA);
         xg = (m1 * lab + m2 * lac * cosA) * invm;
         yg = m2 * lac * sinA * invm;


         // vectors GA, GB, GC
         pos_prec xga, yga, xgb, ygb, xgc, ygc;
         pos_prec lga, lgb, lgc;
         pos_prec cosAGB, cosAGC, sinAGB, sinAGC;
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
      pos_prec sinphi, cosphi, sinpsi, cospsi;
      sinphi = za1d / yra;
      cosphi = sqrt(1 - sinphi * sinphi);
      sinpsi = ((zb1d - zc1d) - (yrb - yrc) * sinphi) / ((xrc - xrb) * cosphi);
      cospsi = sqrt(1 - sinpsi * sinpsi);


      //====================================================================//
      // Eq.A3 local frame vectors da2', db2', dc2'
      pos_prec ya2d;
      pos_prec xb2d, yb2d;
      pos_prec xc2d, yc2d;
      ya2d = yra * cosphi;
      xb2d = xrb * cospsi;
      yb2d = yrb * cosphi + xrb * sinphi * sinpsi;
      xc2d = xrc * cospsi;
      yc2d = yrc * cosphi + xrc * sinphi * sinpsi;
      // adjust numerical error for xb2' and xc2'
      // deltx**2 + ybc**2 + zbc**2 = lbc**2, where ybc**2 + zbc**2 = hh2
      pos_prec deltx, hh2;
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
      pos_prec alpa, beta, gama, al2be2, sinthe, costhe;
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
      pos_prec xa3d, ya3d, za3d;
      pos_prec xb3d, yb3d, zb3d;
      pos_prec xc3d, yc3d, zc3d;
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
      pos_prec xa3, ya3, za3;
      pos_prec xb3, yb3, zb3;
      pos_prec xc3, yc3, zc3;
      xa3 = trns11 * xa3d + trns12 * ya3d + trns13 * za3d;
      ya3 = trns21 * xa3d + trns22 * ya3d + trns23 * za3d;
      za3 = trns31 * xa3d + trns32 * ya3d + trns33 * za3d;
      xb3 = trns11 * xb3d + trns12 * yb3d + trns13 * zb3d;
      yb3 = trns21 * xb3d + trns22 * yb3d + trns23 * zb3d;
      zb3 = trns31 * xb3d + trns32 * yb3d + trns33 * zb3d;
      xc3 = trns11 * xc3d + trns12 * yc3d + trns13 * zc3d;
      yc3 = trns21 * xc3d + trns22 * yc3d + trns23 * zc3d;
      zc3 = trns31 * xc3d + trns32 * yc3d + trns33 * zc3d;


      xpos[ia] = xcom + xa3;
      ypos[ia] = ycom + ya3;
      zpos[ia] = zcom + za3;
      xpos[ib] = xcom + xb3;
      ypos[ib] = ycom + yb3;
      zpos[ib] = zcom + zb3;
      xpos[ic] = xcom + xc3;
      ypos[ic] = ycom + yc3;
      zpos[ic] = zcom + zc3;


      //====================================================================//
      // velocity correction due to the position constraints


      // This code is not in the original SETTLE paper, but is necessary for
      // the velocity verlet integrator.
      vx[ia] += (xa3 - xa1) * invdt;
      vy[ia] += (ya3 - ya1) * invdt;
      vz[ia] += (za3 - za1) * invdt;
      vx[ib] += (xb3 - xb1) * invdt;
      vy[ib] += (yb3 - yb1) * invdt;
      vz[ib] += (zb3 - zb1) * invdt;
      vx[ic] += (xc3 - xc1) * invdt;
      vy[ic] += (yc3 - yc1) * invdt;
      vz[ic] += (zc3 - zc1) * invdt;
   }
}


void rattle2_settle_acc()
{
   if (nratwt <= 0)
      return;


   #pragma acc parallel loop independent async\
           deviceptr(vx,vy,vz,xpos,ypos,zpos,mass,iratwt)
   for (int iw = 0; iw < nratwt; ++iw) {
      int ia, ib, ic;
      mass_prec m0, m1, m2;
      ia = iratwt[iw][0];
      ib = iratwt[iw][1];
      ic = iratwt[iw][2];
      m0 = mass[ia];
      m1 = mass[ib];
      m2 = mass[ic];


      // vectors AB, BC, CA
      pos_prec xab, yab, zab;
      pos_prec xbc, ybc, zbc;
      pos_prec xca, yca, zca;
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
      pos_prec ablng, bclng, calng;
      pos_prec xeab, yeab, zeab;
      pos_prec xebc, yebc, zebc;
      pos_prec xeca, yeca, zeca;
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
      vel_prec xvab, yvab, zvab;
      vel_prec xvbc, yvbc, zvbc;
      vel_prec xvca, yvca, zvca;
      xvab = vx[ib] - vx[ia];
      yvab = vy[ib] - vy[ia];
      zvab = vz[ib] - vz[ia];
      xvbc = vx[ic] - vx[ib];
      yvbc = vy[ic] - vy[ib];
      zvbc = vz[ic] - vz[ib];
      xvca = vx[ia] - vx[ic];
      yvca = vy[ia] - vy[ic];
      zvca = vz[ia] - vz[ic];


      vel_prec vabab, vbcbc, vcaca;
      vabab = xvab * xeab + yvab * yeab + zvab * zeab;
      vbcbc = xvbc * xebc + yvbc * yebc + zvbc * zebc;
      vcaca = xvca * xeca + yvca * yeca + zvca * zeca;


      pos_prec cosa, cosb, cosc;
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
      pos_prec a1, a2, a3, b1, b2, b3, c1, c2, c3, denom;
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
      pos_prec av1, av2, av3, bv1, bv2, bv3, cv1, cv2, cv3;
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
      pos_prec tabd, tbcd, tcad;
      tabd =
         av1 * m0 * m1 * vabab + av2 * m1 * m2 * vbcbc + av3 * m2 * m0 * vcaca;
      tbcd =
         bv1 * m0 * m1 * vabab + bv2 * m1 * m2 * vbcbc + bv3 * m2 * m0 * vcaca;
      tcad =
         cv1 * m0 * m1 * vabab + cv2 * m1 * m2 * vbcbc + cv3 * m2 * m0 * vcaca;


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
   }
}
}
