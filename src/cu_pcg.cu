#include "cudalib.h"
#include "e_polar.h"
#include "induce.h"
#include "io_print.h"
#include "launch.h"
#include "tinker_rt.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/polpcg.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/units.hh>


TINKER_NAMESPACE_BEGIN
#define ITHREAD threadIdx.x + blockIdx.x* blockDim.x
#define STRIDE blockDim.x* gridDim.x


extern "C" __global__
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


extern "C" __global__
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


extern "C" __global__
void pcg_p2(int n, const real* restrict ka, const real* restrict kap,
            const real* restrict ksum, const real* restrict ksump,
            real (*restrict uind)[3], real (*restrict uinp)[3],
            const real (*restrict conj)[3], const real (*restrict conjp)[3],
            real (*restrict rsd)[3], real (*restrict rsdp)[3],
            const real (*restrict vec)[3], const real (*restrict vecp)[3])
{
   real a = *ksum / *ka;
   real ap = *ksump / *kap;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         uind[i][j] += a * conj[i][j];
         uinp[i][j] += ap * conjp[i][j];
         rsd[i][j] -= a * vec[i][j];
         rsdp[i][j] -= ap * vecp[i][j];
      }
   }
}


extern "C" __global__
void pcg_p3(int n, const real* restrict ksum, const real* restrict ksump,
            const real* restrict ksum1, const real* restrict ksump1,
            real (*restrict conj)[3], real (*restrict conjp)[3],
            real (*restrict zrsd)[3], real (*restrict zrsdp)[3])
{
   real b = *ksum1 / *ksum;
   real bp = *ksump1 / *ksump;
   for (int i = ITHREAD; i < n; i += STRIDE) {
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         conj[i][j] = zrsd[i][j] + b * conj[i][j];
         conjp[i][j] = zrsdp[i][j] + bp * conjp[i][j];
      }
   }
}


extern "C" __global__
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


   const bool dirguess = polpcg::pcgguess;
   const bool sparse_prec = polpcg::pcgprec;


   // zero out the induced dipoles at each site
   device_array::zero(false, n, uind, uinp);


   // get the electrostatic field due to permanent multipoles
   dfield(field, fieldp);


   // direct induced dipoles
   launch_k1s(nonblk, n, pcg_udir, n, polarity, udir, udirp, field, fieldp);
   if (dirguess) {
      device_array::copy(false, n, uind, udir);
      device_array::copy(false, n, uinp, udirp);
   }


   // initial residual r(0)
   // if do not use pcgguess, r(0) = E - T Zero = E
   // if use pcgguess, r(0) = E - (inv_alpha + Tu) alpha E
   //                       = E - E -Tu udir
   //                       = -Tu udir
   if (dirguess) {
      ufield(udir, udirp, rsd, rsdp);
   } else {
      device_array::copy(false, n, rsd, field);
      device_array::copy(false, n, rsdp, fieldp);
   }


   // initial M r(0) and p(0)
   if (sparse_prec) {
      sparse_precond_build();
      sparse_precond_apply(rsd, rsdp, zrsd, zrsdp);
   } else {
      diag_precond(rsd, rsdp, zrsd, zrsdp);
   }
   device_array::copy(false, n, conj, zrsd);
   device_array::copy(false, n, conjp, zrsdp);


   // initial r(0) M r(0)
   real* sum = &dptr_real64[0];
   real* sump = &dptr_real64[1];
   device_array::dot(false, n, sum, rsd, zrsd);
   device_array::dot(false, n, sump, rsdp, zrsdp);


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
      real* a = &dptr_real64[2];
      real* ap = &dptr_real64[3];
      // a <- r M r / p T p; a = sum / a; ap = sump / ap
      device_array::dot(false, n, a, conj, vec);
      device_array::dot(false, n, ap, conjp, vecp);


      // u <- u + a p
      // r <- r - a T p
      launch_k1s(nonblk, n, pcg_p2, n, a, ap, sum, sump, uind, uinp, conj,
                 conjp, rsd, rsdp, vec, vecp);


      // calculate/update M r
      if (sparse_prec)
         sparse_precond_apply(rsd, rsdp, zrsd, zrsdp);
      else
         diag_precond(rsd, rsdp, zrsd, zrsdp);


      // b = sum1 / sum; bp = sump1 / sump
      real* sum1 = &dptr_real64[4];
      real* sump1 = &dptr_real64[5];
      device_array::dot(false, n, sum1, rsd, zrsd);
      device_array::dot(false, n, sump1, rsdp, zrsdp);


      // calculate/update p
      launch_k1s(nonblk, n, pcg_p3, n, sum, sump, sum1, sump1, conj, conjp,
                 zrsd, zrsdp);


      // copy sum1/p to sum/p
      device_array::copy(false, 2, sum, sum1);


      real* epsd = &dptr_real64[6];
      real* epsp = &dptr_real64[7];
      device_array::dot(false, n, epsd, rsd, rsd);
      device_array::dot(false, n, epsp, rsdp, rsdp);
      check_rt(cudaMemcpyAsync(pinned_real64, epsd, 2 * sizeof(real),
                               cudaMemcpyDeviceToHost, nonblk));
      check_rt(cudaStreamSynchronize(nonblk));
      epsold = eps;
      eps = REAL_MAX(pinned_real64[0], pinned_real64[1]);
      eps = debye * REAL_SQRT(eps / n);


      if (debug) {
         if (iter == 1) {
            print(stdout,
                  "\n Determination of SCF Induced Dipole Moments\n\n"
                  "{0:4s}Iter{0:4s}RMS Residual (Debye)\n\n",
                  "");
         }
         print(stdout, "{0:>8d}{2:8s}{1:<16.10f}\n", iter, eps, "");
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
            " Induced Dipoles :{2:4s}Iterations{0:>5d}{2:6s}RMS "
            "Residual{1:>15.10f}\n",
            iter, eps, "");
   }


   // terminate the calculation if dipoles failed to converge
   if (iter >= maxiter || eps > epsold) {
      TINKER_RT(prterr)();
      TINKER_THROW("INDUCE  --  Warning, Induced Dipoles are not Converged");
   }
}
TINKER_NAMESPACE_END
