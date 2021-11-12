#include "epolar_chgpen.h"
#include "field_chgpen.h"
#include "induce_donly.h"
#include "launch.h"
#include "mod.uprior.h"
#include "tinker_rt.h"
#include "tool/cudalib.h"
#include "tool/io_print.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/polpcg.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/units.hh>


namespace tinker {
#define ITHREAD threadIdx.x + blockIdx.x* blockDim.x
#define STRIDE  blockDim.x* gridDim.x


__global__
void pcg_udir_donly(int n, const real* restrict polarity,
                    real (*restrict udir)[3], const real (*restrict field)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli = polarity[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         udir[i][j] = poli * field[i][j];
      }
   }
}


__global__
void pcg_rsd2(int n, const real* restrict polarity_inv, //
              real (*restrict rsd)[3],                  //
              const real (*restrict udir)[3], const real (*restrict uind)[3],
              const real (*restrict field)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli_inv = polarity_inv[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j)
         rsd[i][j] = (udir[i][j] - uind[i][j]) * poli_inv + field[i][j];
   }
}


__global__
void pcg_rsd1(int n, const real* restrict polarity, real (*restrict rsd)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      if (polarity[i] == 0) {
         rsd[i][0] = 0;
         rsd[i][1] = 0;
         rsd[i][2] = 0;
      }
   }
}


__global__
void pcg_p4(int n, const real* restrict polarity_inv, real (*restrict vec)[3],
            const real (*restrict conj)[3], const real (*restrict field)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli_inv = polarity_inv[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j)
         vec[i][j] = poli_inv * conj[i][j] - field[i][j];
   }
}


__global__
void pcg_p5(int n, const real* restrict polarity, //
            const real* restrict ka,              //
            const real* restrict ksum, real (*restrict uind)[3],
            const real (*restrict conj)[3], real (*restrict rsd)[3],
            const real (*restrict vec)[3])
{
   real kaval = *ka;
   real a = *ksum / kaval;
   if (kaval == 0)
      a = 0;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      #pragma unroll
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
}


__global__
void pcg_p6(int n, const real* restrict ksum, const real* restrict ksum1,
            real (*restrict conj)[3], real (*restrict zrsd)[3])
{
   real ksumval = *ksum;
   real b = *ksum1 / ksumval;
   if (ksumval == 0)
      b = 0;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      #pragma unroll
      for (int j = 0; j < 3; ++j)
         conj[i][j] = zrsd[i][j] + b * conj[i][j];
   }
}


__global__
void pcg_peek1(int n, float pcgpeek, const real* restrict polarity,
               real (*restrict uind)[3], const real (*restrict rsd)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real term = pcgpeek * polarity[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j)
         uind[i][j] += term * rsd[i][j];
   }
}


void induce_mutual_pcg2_cu(real (*uind)[3])
{
   auto* field = work01_;
   auto* rsd = work02_;
   auto* zrsd = work03_;
   auto* conj = work04_;
   auto* vec = work05_;


   const bool sparse_prec = polpcg::pcgprec;
   bool dirguess = polpcg::pcgguess;
   bool predict = polpred != UPred::NONE;
   if (predict and nualt < maxualt) {
      predict = false;
      dirguess = true;
   }


   // get the electrostatic field due to permanent multipoles
   dfield_chgpen(field);
   // direct induced dipoles
   launch_k1s(g::s0, n, pcg_udir_donly, n, polarity, udir, field);


   // initial induced dipole
   if (predict) {
      ulspred_sum2(uind);
   } else if (dirguess) {
      darray::copy(g::q0, n, uind, udir);
   } else {
      darray::zero(g::q0, n, uind);
   }


   // initial residual r(0)
   //
   // if use pcgguess, r(0) = E - (inv_alpha + Tu) alpha E
   //                       = E - E -Tu udir
   //                       = -Tu udir
   //
   // in general, r(0) = E - (inv_alpha + Tu) u(0)
   //                  = -Tu u(0) + E - inv_alpha u(0)
   //                  = -Tu u(0) + inv_alpha (udir - u(0))
   //
   // if do not use pcgguess, r(0) = E - T Zero = E
   if (predict) {
      ufield_chgpen(uind, field);
      launch_k1s(g::s0, n, pcg_rsd2, n, polarity_inv, rsd, udir, uind, field);
   } else if (dirguess) {
      ufield_chgpen(udir, rsd);
   } else {
      darray::copy(g::q0, n, rsd, field);
   }
   launch_k1s(g::s0, n, pcg_rsd1, n, polarity, rsd);


   // initial M r(0) and p(0)
   if (sparse_prec) {
      sparse_precond_build2();
      sparse_precond_apply2(rsd, zrsd);
   } else {
      diag_precond2(rsd, zrsd);
   }
   darray::copy(g::q0, n, conj, zrsd);


   // initial r(0) M r(0)
   real* sum = &((real*)dptr_buf)[0];
   darray::dot(g::q0, n, sum, rsd, zrsd);


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


      // T p and p
      // vec = (inv_alpha + Tu) conj, field = -Tu conj
      // vec = inv_alpha * conj - field
      ufield_chgpen(conj, field);
      launch_k1s(g::s0, n, pcg_p4, n, polarity_inv, vec, conj, field);


      // a <- p T p
      real* a = &((real*)dptr_buf)[1];
      // a <- r M r / p T p; a = sum / a; ap = sump / ap
      darray::dot(g::q0, n, a, conj, vec);


      // u <- u + a p
      // r <- r - a T p
      launch_k1s(g::s0, n, pcg_p5, n, polarity, a, sum, uind, conj, rsd, vec);


      // calculate/update M r
      if (sparse_prec)
         sparse_precond_apply2(rsd, zrsd);
      else
         diag_precond2(rsd, zrsd);


      // b = sum1 / sum; bp = sump1 / sump
      real* sum1 = &((real*)dptr_buf)[2];
      darray::dot(g::q0, n, sum1, rsd, zrsd);


      // calculate/update p
      launch_k1s(g::s0, n, pcg_p6, n, sum, sum1, conj, zrsd);


      // copy sum1/p to sum/p
      darray::copy(g::q0, 2, sum, sum1);


      real* epsd = &((real*)dptr_buf)[3];
      darray::dot(g::q0, n, epsd, rsd, rsd);
      check_rt(cudaMemcpyAsync((real*)pinned_buf, epsd, sizeof(real),
                               cudaMemcpyDeviceToHost, g::s0));
      check_rt(cudaStreamSynchronize(g::s0));
      epsold = eps;
      eps = ((real*)pinned_buf)[0];
      eps = debye * REAL_SQRT(eps / n);


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
      if (done)
         launch_k1s(g::s0, n, pcg_peek1, n, pcgpeek, polarity, uind, rsd);
   }


   // print the results from the conjugate gradient iteration
   if (debug) {
      print(stdout,
            " Induced Dipoles :    Iterations %4d      RMS"
            " Residual %14.10f\n",
            iter, eps);
   }


   // terminate the calculation if dipoles failed to converge
   if (iter >= maxiter || eps > epsold) {
      prterr();
      TINKER_THROW("INDUCE  --  Warning, Induced Dipoles are not Converged");
   }
}
}
