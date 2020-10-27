#include "epolar.h"
#include "field.h"
#include "induce.h"
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
void pcg_udir(int n, const real* restrict polarity, real (*restrict udir)[3],
              real (*restrict udirp)[3], const real (*restrict field)[3],
              const real (*restrict fieldp)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli = polarity[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         udir[i][j] = poli * field[i][j];
         udirp[i][j] = poli * fieldp[i][j];
      }
   }
}


__global__
void pcg_rsd(int n, const real* restrict polarity_inv,         //
             real (*restrict rsd)[3], real (*restrict rsp)[3], //
             const real (*restrict udir)[3], const real (*restrict udip)[3],
             const real (*restrict uind)[3], const real (*restrict uinp)[3],
             const real (*restrict field)[3], const real (*restrict fielp)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli_inv = polarity_inv[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         rsd[i][j] = (udir[i][j] - uind[i][j]) * poli_inv + field[i][j];
         rsp[i][j] = (udip[i][j] - uinp[i][j]) * poli_inv + fielp[i][j];
      }
   }
}


__global__
void pcg_rsd0(int n, const real* restrict polarity, real (*restrict rsd)[3],
              real (*restrict rsdp)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      if (polarity[i] == 0) {
         rsd[i][0] = 0;
         rsd[i][1] = 0;
         rsd[i][2] = 0;
         rsdp[i][0] = 0;
         rsdp[i][1] = 0;
         rsdp[i][2] = 0;
      }
   }
}


__global__
void pcg_p1(int n, const real* restrict polarity_inv, real (*restrict vec)[3],
            real (*restrict vecp)[3], const real (*restrict conj)[3],
            const real (*restrict conjp)[3], const real (*restrict field)[3],
            const real (*restrict fieldp)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli_inv = polarity_inv[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         vec[i][j] = poli_inv * conj[i][j] - field[i][j];
         vecp[i][j] = poli_inv * conjp[i][j] - fieldp[i][j];
      }
   }
}


__global__
void pcg_p2(int n, const real* restrict polarity,              //
            const real* restrict ka, const real* restrict kap, //
            const real* restrict ksum, const real* restrict ksump,
            real (*restrict uind)[3], real (*restrict uinp)[3],
            const real (*restrict conj)[3], const real (*restrict conjp)[3],
            real (*restrict rsd)[3], real (*restrict rsdp)[3],
            const real (*restrict vec)[3], const real (*restrict vecp)[3])
{
   real kaval = *ka, kapval = *kap;
   real a = *ksum / kaval, ap = *ksump / kapval;
   if (kaval == 0)
      a = 0;
   if (kapval == 0)
      ap = 0;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         uind[i][j] += a * conj[i][j];
         uinp[i][j] += ap * conjp[i][j];
         rsd[i][j] -= a * vec[i][j];
         rsdp[i][j] -= ap * vecp[i][j];
      }
      if (polarity[i] == 0) {
         rsd[i][0] = 0;
         rsd[i][1] = 0;
         rsd[i][2] = 0;
         rsdp[i][0] = 0;
         rsdp[i][1] = 0;
         rsdp[i][2] = 0;
      }
   }
}


__global__
void pcg_p3(int n, const real* restrict ksum, const real* restrict ksump,
            const real* restrict ksum1, const real* restrict ksump1,
            real (*restrict conj)[3], real (*restrict conjp)[3],
            real (*restrict zrsd)[3], real (*restrict zrsdp)[3])
{
   real kaval = *ksum, kapval = *ksump;
   real b = *ksum1 / kaval, bp = *ksump1 / kapval;
   if (kaval == 0)
      b = 0;
   if (kapval == 0)
      bp = 0;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         conj[i][j] = zrsd[i][j] + b * conj[i][j];
         conjp[i][j] = zrsdp[i][j] + bp * conjp[i][j];
      }
   }
}


__global__
void pcg_peek(int n, float pcgpeek, const real* restrict polarity,
              real (*restrict uind)[3], real (*restrict uinp)[3],
              const real (*restrict rsd)[3], const real (*restrict rsdp)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real term = pcgpeek * polarity[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         uind[i][j] += term * rsd[i][j];
         uinp[i][j] += term * rsdp[i][j];
      }
   }
}


void induce_mutual_pcg1_cu(real (*uind)[3], real (*uinp)[3])
{
   auto* field = work01_;
   auto* fieldp = work02_;
   auto* rsd = work03_;
   auto* rsdp = work04_;
   auto* zrsd = work05_;
   auto* zrsdp = work06_;
   auto* conj = work07_;
   auto* conjp = work08_;
   auto* vec = work09_;
   auto* vecp = work10_;


   const bool sparse_prec = polpcg::pcgprec;
   bool dirguess = polpcg::pcgguess;
   bool predict = polpred != UPred::NONE;
   if (predict and nualt < maxualt) {
      predict = false;
      dirguess = true;
   }


   // get the electrostatic field due to permanent multipoles
   dfield(field, fieldp);
   // direct induced dipoles
   launch_k1s(nonblk, n, pcg_udir, n, polarity, udir, udirp, field, fieldp);


   // initial induced dipole
   if (predict) {
      ulspred_sum(uind, uinp);
   } else if (dirguess) {
      darray::copy(PROCEED_NEW_Q, n, uind, udir);
      darray::copy(PROCEED_NEW_Q, n, uinp, udirp);
   } else {
      darray::zero(PROCEED_NEW_Q, n, uind, uinp);
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
      ufield(uind, uinp, field, fieldp);
      launch_k1s(nonblk, n, pcg_rsd, n, polarity_inv, rsd, rsdp, udir, udirp,
                 uind, uinp, field, fieldp);
   } else if (dirguess) {
      ufield(udir, udirp, rsd, rsdp);
   } else {
      darray::copy(PROCEED_NEW_Q, n, rsd, field);
      darray::copy(PROCEED_NEW_Q, n, rsdp, fieldp);
   }
   launch_k1s(nonblk, n, pcg_rsd0, n, polarity, rsd, rsdp);


   // initial M r(0) and p(0)
   if (sparse_prec) {
      sparse_precond_build();
      sparse_precond_apply(rsd, rsdp, zrsd, zrsdp);
   } else {
      diag_precond(rsd, rsdp, zrsd, zrsdp);
   }
   darray::copy(PROCEED_NEW_Q, n, conj, zrsd);
   darray::copy(PROCEED_NEW_Q, n, conjp, zrsdp);


   // initial r(0) M r(0)
   real* sum = &((real*)dptr_buf)[0];
   real* sump = &((real*)dptr_buf)[1];
   darray::dot(PROCEED_NEW_Q, n, sum, rsd, zrsd);
   darray::dot(PROCEED_NEW_Q, n, sump, rsdp, zrsdp);


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
      ufield(conj, conjp, field, fieldp);
      launch_k1s(nonblk, n, pcg_p1, n, polarity_inv, vec, vecp, conj, conjp,
                 field, fieldp);


      // a <- p T p
      real* a = &((real*)dptr_buf)[2];
      real* ap = &((real*)dptr_buf)[3];
      // a <- r M r / p T p; a = sum / a; ap = sump / ap
      darray::dot(PROCEED_NEW_Q, n, a, conj, vec);
      darray::dot(PROCEED_NEW_Q, n, ap, conjp, vecp);


      // u <- u + a p
      // r <- r - a T p
      launch_k1s(nonblk, n, pcg_p2, n, polarity, a, ap, sum, sump, uind, uinp,
                 conj, conjp, rsd, rsdp, vec, vecp);


      // calculate/update M r
      if (sparse_prec)
         sparse_precond_apply(rsd, rsdp, zrsd, zrsdp);
      else
         diag_precond(rsd, rsdp, zrsd, zrsdp);


      // b = sum1 / sum; bp = sump1 / sump
      real* sum1 = &((real*)dptr_buf)[4];
      real* sump1 = &((real*)dptr_buf)[5];
      darray::dot(PROCEED_NEW_Q, n, sum1, rsd, zrsd);
      darray::dot(PROCEED_NEW_Q, n, sump1, rsdp, zrsdp);


      // calculate/update p
      launch_k1s(nonblk, n, pcg_p3, n, sum, sump, sum1, sump1, conj, conjp,
                 zrsd, zrsdp);


      // copy sum1/p to sum/p
      darray::copy(PROCEED_NEW_Q, 2, sum, sum1);


      real* epsd = &((real*)dptr_buf)[6];
      real* epsp = &((real*)dptr_buf)[7];
      darray::dot(PROCEED_NEW_Q, n, epsd, rsd, rsd);
      darray::dot(PROCEED_NEW_Q, n, epsp, rsdp, rsdp);
      check_rt(cudaMemcpyAsync((real*)pinned_buf, epsd, 2 * sizeof(real),
                               cudaMemcpyDeviceToHost, nonblk));
      check_rt(cudaStreamSynchronize(nonblk));
      epsold = eps;
      eps = REAL_MAX(((real*)pinned_buf)[0], ((real*)pinned_buf)[1]);
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
      if (done) {
         launch_k1s(nonblk, n, pcg_peek, n, pcgpeek, polarity, uind, uinp, rsd,
                    rsdp);
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
}
