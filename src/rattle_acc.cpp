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
}
