#include "add.h"
#include "epolar.h"
#include "empole_chgpen.h"
#include "epolar_chgpen.h"
#include "field_chgpen.h"
#include "glob.nblist.h"
#include "image.h"
#include "induce_donly.h"
#include "seq_damp_chgpen.h"
#include "tinker_rt.h"
#include "tool/error.h"
#include "tool/gpu_card.h"
#include "tool/io_print.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/polpcg.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/units.hh>

namespace tinker {
void diag_precond2(const real (*rsd)[3], real (*zrsd)[3])
{
   #pragma acc parallel loop independent async\
               deviceptr(polarity,rsd,zrsd)
   for (int i = 0; i < n; ++i) {
      real poli = polarity[i];
      #pragma acc loop seq
      for (int j = 0; j < 3; ++j) 
         zrsd[i][j] = poli * rsd[i][j];
   }
}

#define APPLY_DPTRS rsd, zrsd, x, y, z, polarity, palpha
void sparse_precond_apply_acc2(const real (*rsd)[3], 
                              real (*zrsd)[3])
{
   #pragma acc parallel loop independent async\
               deviceptr(polarity,rsd,zrsd)
   for (int i = 0; i < n; ++i) {
      real poli = udiag * polarity[i];
      #pragma acc loop seq
      for (int j = 0; j < 3; ++j) 
         zrsd[i][j] = poli * rsd[i][j];
   }

   const int maxnlst = ulist_unit->maxnlst;
   const auto* ulst = ulist_unit.deviceptr();

   MAYBE_UNUSED int GRID_DIM = get_grid_size(BLOCK_DIM);
   #pragma acc parallel async num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               deviceptr(APPLY_DPTRS,ulst)
   #pragma acc loop gang independent
   for (int i = 0; i < n; ++i) {
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real alphai = palpha[i];
      real poli = polarity[i];

      int nulsti = ulst->nlst[i];
      int base = i * maxnlst;
      real gxi = 0, gyi = 0, gzi = 0;

      #pragma acc loop vector independent reduction(+:gxi,gyi,gzi)
      for (int kk = 0; kk < nulsti; ++kk) {
         int k = ulst->lst[base + kk];

         real alphak = palpha[k];
         real xr = x[k] - xi;
         real yr = y[k] - yi;
         real zr = z[k] - zi;

         real r2 = image2(xr, yr, zr);
         real r = REAL_SQRT(r2);

         real dmpik[3];
         damp_mut(dmpik,r,alphai,alphak);
         real scale3 = dmpik[1];
         real scale5 = dmpik[2];

         real polik = poli * polarity[k];
         real rr3 = scale3 * polik * REAL_RECIP(r * r2);
         real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);

         real m0 = rr5 * xr * xr - rr3;
         real m1 = rr5 * xr * yr;
         real m2 = rr5 * xr * zr;
         real m3 = rr5 * yr * yr - rr3;
         real m4 = rr5 * yr * zr;
         real m5 = rr5 * zr * zr - rr3;

         gxi += m0 * rsd[k][0] + m1 * rsd[k][1] + m2 * rsd[k][2];
         gyi += m1 * rsd[k][0] + m3 * rsd[k][1] + m4 * rsd[k][2];
         gzi += m2 * rsd[k][0] + m4 * rsd[k][1] + m5 * rsd[k][2];

         atomic_add(m0 * rsd[i][0] + m1 * rsd[i][1] + m2 * rsd[i][2],
                    &zrsd[k][0]);
         atomic_add(m1 * rsd[i][0] + m3 * rsd[i][1] + m4 * rsd[i][2],
                    &zrsd[k][1]);
         atomic_add(m2 * rsd[i][0] + m4 * rsd[i][1] + m5 * rsd[i][2],
                    &zrsd[k][2]);

      }

      atomic_add(gxi, &zrsd[i][0]);
      atomic_add(gyi, &zrsd[i][1]);
      atomic_add(gzi, &zrsd[i][2]);
   }

   #pragma acc parallel loop async\
               deviceptr(APPLY_DPTRS,wexclude,wexclude_scale)
   for (int ii = 0; ii < nwexclude; ++ii) {
      int i = wexclude[ii][0];
      int k = wexclude[ii][1];
      real wscale = wexclude_scale[ii];

      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real alphai = palpha[i];
      real alphak = palpha[k];

      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;

      real r2 = image2(xr, yr, zr);
      real r = REAL_SQRT(r2);

      real dmpik[3];
      damp_mut(dmpik,r,alphai,alphak);
      real scale3 = wscale * dmpik[1];
      real scale5 = wscale * dmpik[2];

      real poli = polarity[i];
      real polik = poli * polarity[k];
      real rr3 = scale3 * polik * REAL_RECIP(r * r2);
      real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);

      real m0 = rr5 * xr * xr - rr3;
      real m1 = rr5 * xr * yr;
      real m2 = rr5 * xr * zr;
      real m3 = rr5 * yr * yr - rr3;
      real m4 = rr5 * yr * zr;
      real m5 = rr5 * zr * zr - rr3;

      atomic_add(m0 * rsd[k][0] + m1 * rsd[k][1] + m2 * rsd[k][2], &zrsd[i][0]);
      atomic_add(m1 * rsd[k][0] + m3 * rsd[k][1] + m4 * rsd[k][2], &zrsd[i][1]);
      atomic_add(m2 * rsd[k][0] + m4 * rsd[k][1] + m5 * rsd[k][2], &zrsd[i][2]);

      atomic_add(m0 * rsd[i][0] + m1 * rsd[i][1] + m2 * rsd[i][2], &zrsd[k][0]);
      atomic_add(m1 * rsd[i][0] + m3 * rsd[i][1] + m4 * rsd[i][2], &zrsd[k][1]);
      atomic_add(m2 * rsd[i][0] + m4 * rsd[i][1] + m5 * rsd[i][2], &zrsd[k][2]);
   }
}

/**
 * PCG
 *
 * M = preconditioner
 * T = inv_alpha + Tu, -Tu = ufield_chgpen
 *
 * r(0) = E - T u(0)
 * p(0) (or conj(0)) = M r(0)
 * ------------------------------
 * gamma (or a) = r(i) M r(i) / p(i) T p(i)
 * u(i+1) = u(i) + a p(i)
 * r(i+1) = r(i) - a T p(i)
 * beta (or b(i+1)) = r(i+1) M r(i+1) / r(i) M r(i)
 * p(i+1) = M r(r+1) + b(i+1) p(i)
 * ------------------------------
 *
 * subroutine induce0a in induce.f
 * rsd = r
 * zrsd = M r
 * conj = p
 * vec = T P
 */
void induce_mutual_pcg_acc2(real (*uind)[3])
{
   auto* field = work01_;
   auto* rsd = work02_;
   auto* zrsd = work03_;
   auto* conj = work04_;
   auto* vec = work05_;

   const bool dirguess = polpcg::pcgguess;
   // use sparse matrix preconditioner
   // or just use diagonal matrix preconditioner
   const bool sparse_prec = polpcg::pcgprec;

   // zero out the induced dipoles at each site

   darray::zero(PROCEED_NEW_Q, n, uind);

   // get the electrostatic field due to permanent multipoles

   dfield_chgpen(field);

   // direct induced dipoles

   #pragma acc parallel loop independent async\
               deviceptr(polarity,udir,field)
   for (int i = 0; i < n; ++i) {
      real poli = polarity[i];
      #pragma acc loop seq
      for (int j = 0; j < 3; ++j) 
         udir[i][j] = poli * field[i][j];
   }

   if (dirguess) 
      darray::copy(PROCEED_NEW_Q, n, uind, udir);

   // initial residual r(0)

   // if do not use pcgguess, r(0) = E - T Zero = E
   // if use pcgguess, r(0) = E - (inv_alpha + Tu) alpha E
   //                       = E - E -Tu udir
   //                       = -Tu udir

   if (dirguess) {
      ufield_chgpen(udir, rsd);
   } else {
      darray::copy(PROCEED_NEW_Q, n, rsd, field);
   }
   #pragma acc parallel loop independent async deviceptr(polarity,rsd)
   for (int i = 0; i < n; ++i) {
      if (polarity[i] == 0) {
         rsd[i][0] = 0;
         rsd[i][1] = 0;
         rsd[i][2] = 0;
      }
   }

   // initial M r(0) and p(0)

   if (sparse_prec) {
      sparse_precond_build2();
      sparse_precond_apply2(rsd, zrsd);
   } else {
      diag_precond2(rsd, zrsd);
   }
   darray::copy(PROCEED_NEW_Q, n, conj, zrsd);

   // initial r(0) M r(0)

   real sum;
   sum = darray::dot(WAIT_NEW_Q, n, rsd, zrsd);

   // conjugate gradient iteration of the mutual induced dipoles

   const bool debug = inform::debug;
   const int politer = polpot::politer;
   const real poleps = polpot::poleps;
   const real debye = units::debye;
   const real pcgpeek = polpcg::pcgpeek;
   const int maxiter = 100; // see also subroutine induce0a in induce.f

   bool done = false;
   int iter = 0;
   real eps = 100;
   real epsold;

   while (!done) {
      ++iter;

      // vec = (inv_alpha + Tu) conj, field = -Tu conj
      // vec = inv_alpha * conj - field
      ufield_chgpen(conj, field);
      #pragma acc parallel loop independent async\
                  deviceptr(polarity_inv,vec,conj,field)
      for (int i = 0; i < n; ++i) {
         real poli_inv = polarity_inv[i];
         #pragma acc loop seq
         for (int j = 0; j < 3; ++j) 
            vec[i][j] = poli_inv * conj[i][j] - field[i][j];
      }

      // a <- p T p
      real a;
      a = darray::dot(WAIT_NEW_Q, n, conj, vec);
      // a <- r M r / p T p
      if (a != 0)
         a = sum / a;

      // u <- u + a p
      // r <- r - a T p
      #pragma acc parallel loop independent async\
                  deviceptr(polarity,uind,conj,rsd,vec)
      for (int i = 0; i < n; ++i) {
         #pragma acc loop seq
         for (int j = 0; j < 3; ++j) {
            uind[i][j] += a * conj[i][j];
            rsd[i][j] -= a * vec[i][j];
         }
         if (polarity[i] == 0) {
            rsd[i][0] = 0;
            rsd[i][1] = 0;
            rsd[i][2] = 0;
         }
      }

      // calculate/update M r
      if (sparse_prec)
         sparse_precond_apply2(rsd, zrsd);
      else
         diag_precond2(rsd, zrsd);

      real b;
      real sum1;
      sum1 = darray::dot(WAIT_NEW_Q, n, rsd, zrsd);
      if (sum != 0)
         b = sum1 / sum;


      // calculate/update p
      #pragma acc parallel loop independent async\
                  deviceptr(conj,zrsd)
      for (int i = 0; i < n; ++i) {
         #pragma acc loop seq
         for (int j = 0; j < 3; ++j) 
            conj[i][j] = zrsd[i][j] + b * conj[i][j];
      }


      sum = sum1;

      real epsd;
      epsd = darray::dot(WAIT_NEW_Q, n, rsd, rsd);

      epsold = eps;
      eps = epsd;

      if (debug) {
         if (iter == 1) {
            print(stdout,
                  "\n Determination of SCF Induced Dipole Moments\n\n"
                  "    Iter    RMS Residual (Debye)\n\n");
         }
         print(stdout, " %8d       %-16.10f\n", iter, eps);
      }

      if (eps < poleps)
         done = true;
      if (eps > epsold)
         done = true;
      if (iter >= politer)
         done = true;

      // apply a "peek" iteration to the mutual induced dipoles

      if (done) {
         #pragma acc parallel loop independent async\
                     deviceptr(polarity,uind,rsd)
         for (int i = 0; i < n; ++i) {
            real term = pcgpeek * polarity[i];
            #pragma acc loop seq
            for (int j = 0; j < 3; ++j) 
               uind[i][j] += term * rsd[i][j];
         }
      }
   }

   // print the results from the conjugate gradient iteration

   if (debug) {
      print(stdout,
            " Induced Dipoles :    Iterations %4d      RMS "
            "Residual %14.10f\n",
            iter, eps);
   }

   // terminate the calculation if dipoles failed to converge

   if (iter >= maxiter || eps > epsold) {
      t_prterr();
      TINKER_THROW("INDUCE  --  Warning, Induced Dipoles are not Converged");
   }
}


void induce_mutual_pcg2(real (*uind)[3])
{
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      induce_mutual_pcg_cu2(uind);
   else
#endif
      induce_mutual_pcg_acc2(uind);
}
}
