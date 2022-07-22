#include "ff/amoeba/induce.h"
#include "ff/amoebamod.h"
#include "ff/atom.h"
#include "ff/hippo/induce.h"
#include "ff/hippomod.h"
#include "ff/switch.h"
#include "tool/darray.h"
#include "tool/error.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/polpcg.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/units.hh>

namespace tinker {
void induceMutualPcg4_acc(real (*uind)[3])
{
   auto* field = work01_;
   auto* rsd = work02_;
   auto* zrsd = work03_;
   auto* conj = work04_;
   auto* vec = work05_;

   bool dirguess = polpcg::pcgguess;
   // use sparse matrix preconditioner
   // or just use diagonal matrix preconditioner
   const bool sparse_prec = polpcg::pcgprec and (switchOff(Switch::USOLVE) > 0);

   // zero out the induced dipoles at each site

   bool predict = polpred != UPred::NONE;
   if (predict and nualt < maxualt) {
      predict = false;
      dirguess = true;
   }
   // get the electrostatic field due to permanent multipoles
   dfieldChgpen(field);
   // direct induced dipoles

   #pragma acc parallel loop independent async\
               deviceptr(polarity,udir,field)
   for (int i = 0; i < n; ++i) {
      real poli = polarity[i];
      #pragma acc loop seq
      for (int j = 0; j < 3; ++j)
         udir[i][j] = poli * field[i][j];
   }

   alterpol(polscale, polinv);

   // initial induced dipole
   if (predict) {
      ulspredSum(uind, nullptr);
   } else if (dirguess) {
      #pragma acc parallel loop independent async\
              deviceptr(polarity,field,polinv,uind)
      for (int i = 0; i < n; ++i) {
         real poli = polarity[i];
         #pragma acc loop seq
         for (int j = 0; j < 3; ++j) {
            uind[i][j] = poli *
               (polinv[i][0][j] * field[i][0] + polinv[i][1][j] * field[i][1] +
                  polinv[i][2][j] * field[i][2]);
         }
      }
   } else {
      darray::zero(g::q0, n, uind);
   }
   // initial residual r(0)

   // if do not use pcgguess, r(0) = E - T Zero = E
   // if use pcgguess, r(0) = E - (inv_alpha + Tu) alpha E
   //                       = E - E -Tu udir
   //                       = -Tu udir

   if (predict) {
      ufieldChgpen(uind, field);
      #pragma acc parallel loop independent async\
              deviceptr(polarity_inv,udir,uind,field,polscale,rsd)
      for (int i = 0; i < n; ++i) {
         real pol = polarity_inv[i];
         #pragma acc loop seq
         for (int j = 0; j < 3; ++j) {
            rsd[i][j] = (udir[i][j] - uind[i][0] * polscale[i][0][j] -
                           uind[i][1] * polscale[i][1][j] - uind[i][2] * polscale[i][2][j]) *
                  pol +
               field[i][j];
         }
      }
   } else if (dirguess) {
      // uind is used here instead of udir since without exchange polarization udir = uind
      // but with exchange polarization udir != uind (for dirguess).
      ufieldChgpen(uind, rsd);
   } else {
      darray::copy(g::q0, n, rsd, field);
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
      sparsePrecondBuild2();
      sparsePrecondApply2(rsd, zrsd);
   } else {
      diagPrecond2(rsd, zrsd);
   }
   darray::copy(g::q0, n, conj, zrsd);

   // initial r(0) M r(0)

   real sum;
   sum = darray::dotThenReturn(g::q0, n, rsd, zrsd);

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

   while (not done) {
      ++iter;

      // vec = (inv_alpha + Tu) conj, field = -Tu conj
      // vec = inv_alpha * conj - field
      ufieldChgpen(conj, field);
      #pragma acc parallel loop independent async\
                  deviceptr(polarity_inv,conj,field,polscale,vec)
      for (int i = 0; i < n; ++i) {
         real poli_inv = polarity_inv[i];
         #pragma acc loop seq
         for (int j = 0; j < 3; ++j)
            vec[i][j] = poli_inv *
                  (conj[i][0] * polscale[i][0][j] + conj[i][1] * polscale[i][1][j] +
                     conj[i][2] * polscale[i][2][j]) -
               field[i][j];
      }

      // a <- p T p
      real a;
      a = darray::dotThenReturn(g::q0, n, conj, vec);
      // a <- r M r / p T p
      if (a != 0)
         a = sum / a;

      // u <- u + a p
      // r <- r - a T p
      #pragma acc parallel loop independent async\
                  deviceptr(polarity,conj,vec,uind,rsd)
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
         sparsePrecondApply2(rsd, zrsd);
      else
         diagPrecond2(rsd, zrsd);

      real b;
      real sum1;
      sum1 = darray::dotThenReturn(g::q0, n, rsd, zrsd);
      b = sum1 / sum;
      if (sum == 0)
         b = 0;

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
      epsd = darray::dotThenReturn(g::q0, n, rsd, rsd);
      epsold = eps;
      eps = epsd;
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
      // if (eps > epsold)
      //    done = true;
      if (iter >= politer)
         done = true;

      // apply a "peek" iteration to the mutual induced dipoles

      if (done) {
         #pragma acc parallel loop independent async\
                     deviceptr(polarity,rsd,uind)
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

   // if (iter >= maxiter || eps > epsold) {
   if (iter >= maxiter) {
      printError();
      TINKER_THROW("INDUCE  --  Warning, Induced Dipoles are not Converged");
   }
}
}
