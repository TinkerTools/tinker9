#include "ff/amoeba/induce.h"
#include "ff/amoebamod.h"
#include "ff/aplus/induce.h"
#include "ff/atom.h"
#include "ff/image.h"
#include "ff/nblist.h"
#include "ff/switch.h"
#include "seq/add.h"
#include "seq/damp.h"
#include "tool/error.h"
#include "tool/gpucard.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/polpcg.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/units.hh>

namespace tinker {
void induceMutualPcg3_acc(real (*uind)[3])
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
   dfieldAplus(field);
   // direct induced dipoles

   #pragma acc parallel loop independent async\
               deviceptr(polarity,udir,field)
   for (int i = 0; i < n; ++i) {
      real poli = polarity[i];
      #pragma acc loop seq
      for (int j = 0; j < 3; ++j)
         udir[i][j] = poli * field[i][j];
   }

   // initial induced dipole
   if (predict) {
      ulspredSum(uind, nullptr);
   } else if (dirguess) {
      darray::copy(g::q0, n, uind, udir);
   } else {
      darray::zero(g::q0, n, uind);
   }
   // initial residual r(0)

   // if do not use pcgguess, r(0) = E - T Zero = E
   // if use pcgguess, r(0) = E - (inv_alpha + Tu) alpha E
   //                       = E - E -Tu udir
   //                       = -Tu udir

   if (predict) {
      ufieldAplus(uind, field);
      #pragma acc parallel loop independent async\
              deviceptr(polarity_inv,rsd,udir,uind,field)
      for (int i = 0; i < n; ++i) {
         real pol = polarity_inv[i];
         #pragma acc loop seq
         for (int j = 0; j < 3; ++j)
            rsd[i][j] = (udir[i][j] - uind[i][j]) * pol + field[i][j];
      }
   } else if (dirguess) {
      ufieldAplus(udir, rsd);
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
      sparsePrecondBuild3();
      sparsePrecondApply3(rsd, zrsd);
   } else {
      diagPrecond3(rsd, zrsd);
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
      ufieldAplus(conj, field);
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
      a = darray::dotThenReturn(g::q0, n, conj, vec);
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
         sparsePrecondApply3(rsd, zrsd);
      else
         diagPrecond3(rsd, zrsd);

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
      printError();
      TINKER_THROW("INDUCE  --  Warning, Induced Dipoles are not Converged");
   }
}
}

namespace tinker {
#define APPLY_DPTRS rsd, zrsd, x, y, z, polarity, pdamp, thole
void sparsePrecondApply3_acc(const real (*rsd)[3], real (*zrsd)[3])
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

   MAYBE_UNUSED int GRID_DIM = gpuGridSize(BLOCK_DIM);
   #pragma acc parallel async num_gangs(GRID_DIM) vector_length(BLOCK_DIM)\
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
               deviceptr(APPLY_DPTRS,ulst)
   #pragma acc loop gang independent
   for (int i = 0; i < n; ++i) {
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real pdi = pdamp[i];
      real pti = thole[i];
      real poli = polarity[i];

      int nulsti = ulst->nlst[i];
      int base = i * maxnlst;
      real gxi = 0, gyi = 0, gzi = 0;

      #pragma acc loop vector independent reduction(+:gxi,gyi,gzi)
      for (int kk = 0; kk < nulsti; ++kk) {
         int k = ulst->lst[base + kk];

         real pdk = pdamp[k];
         real ptk = thole[k];
         real xr = x[k] - xi;
         real yr = y[k] - yi;
         real zr = z[k] - zi;

         real r2 = image2(xr, yr, zr);
         real r = REAL_SQRT(r2);

         real scale3, scale5;
         damp_thole2(r, pdi, pti, pdk, ptk, scale3, scale5);

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

         atomic_add(m0 * rsd[i][0] + m1 * rsd[i][1] + m2 * rsd[i][2], &zrsd[k][0]);
         atomic_add(m1 * rsd[i][0] + m3 * rsd[i][1] + m4 * rsd[i][2], &zrsd[k][1]);
         atomic_add(m2 * rsd[i][0] + m4 * rsd[i][1] + m5 * rsd[i][2], &zrsd[k][2]);
      }

      atomic_add(gxi, &zrsd[i][0]);
      atomic_add(gyi, &zrsd[i][1]);
      atomic_add(gzi, &zrsd[i][2]);
   }

   #pragma acc parallel loop async\
               present(lvec1,lvec2,lvec3,recipa,recipb,recipc)\
               deviceptr(APPLY_DPTRS,uexclude,uexclude_scale)
   for (int ii = 0; ii < nuexclude; ++ii) {
      int i = uexclude[ii][0];
      int k = uexclude[ii][1];
      real uscale = uexclude_scale[ii] - 1;

      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real pdi = pdamp[i];
      real pti = thole[i];
      real poli = polarity[i];

      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;

      real r2 = image2(xr, yr, zr);
      real r = REAL_SQRT(r2);

      real scale3, scale5;
      damp_thole2(r, pdi, pti, pdamp[k], thole[k], scale3, scale5);
      scale3 *= uscale;
      scale5 *= uscale;

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
}
