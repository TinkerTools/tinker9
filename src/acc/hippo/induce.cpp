#include "ff/hippo/induce.h"
#include "add.h"
#include "ff/amoeba/elecamoeba.h"
#include "ff/atom.h"
#include "ff/hippo/elechippo.h"
#include "ff/image.h"
#include "ff/nblist.h"
#include "math/lu.h"
#include "seq/damp_hippo.h"
#include "tool/error.h"
#include "tool/gpucard.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/polpcg.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/units.hh>

namespace tinker { // PCG
//
// M = preconditioner
// T = inv_alpha + Tu, -Tu = ufield_chgpen
//
// r(0) = E - T u(0)
// p(0) (or conj(0)) = M r(0)
// ------------------------------
// gamma (or a) = r(i) M r(i) / p(i) T p(i)
// u(i+1) = u(i) + a p(i)
// r(i+1) = r(i) - a T p(i)
// beta (or b(i+1)) = r(i+1) M r(i+1) / r(i) M r(i)
// p(i+1) = M r(r+1) + b(i+1) p(i)
// ------------------------------
//
// subroutine induce0a in induce.f
// rsd = r
// zrsd = M r
// conj = p
// vec = T P
void induceMutualPcg2_acc(real (*uind)[3])
{
   auto* field = work01_;
   auto* rsd = work02_;
   auto* zrsd = work03_;
   auto* conj = work04_;
   auto* vec = work05_;

   bool dirguess = polpcg::pcgguess;
   // use sparse matrix preconditioner
   // or just use diagonal matrix preconditioner
   const bool sparse_prec = polpcg::pcgprec;

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

   // initial induced dipole
   if (predict) {
      ulspredSum2(uind);
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
      ufieldChgpen(uind, field);
      #pragma acc parallel loop independent async\
              deviceptr(polarity_inv,rsd,udir,uind,field)
      for (int i = 0; i < n; ++i) {
         real pol = polarity_inv[i];
         #pragma acc loop seq
         for (int j = 0; j < 3; ++j)
            rsd[i][j] = (udir[i][j] - uind[i][j]) * pol + field[i][j];
      }
   } else if (dirguess) {
      ufieldChgpen(udir, rsd);
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
void diagPrecond2_acc(const real (*rsd)[3], real (*zrsd)[3])
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
void sparsePrecondApply2_acc(const real (*rsd)[3], real (*zrsd)[3])
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
         damp_mut(dmpik, r, alphai, alphak);
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
               deviceptr(APPLY_DPTRS,wexclude,wexclude_scale)
   for (int ii = 0; ii < nwexclude; ++ii) {
      int i = wexclude[ii][0];
      int k = wexclude[ii][1];
      real wscale = wexclude_scale[ii] - 1;

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
      damp_mut(dmpik, r, alphai, alphak);
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

void ulspredSave2_acc(const real (*restrict uind)[3])
{
   if (polpred == UPred::NONE)
      return;

   // clang-format off
   real(*ud)[3];
   int pos = nualt % maxualt;
   switch (pos) {
      case  0: ud = udalt_00; break;
      case  1: ud = udalt_01; break;
      case  2: ud = udalt_02; break;
      case  3: ud = udalt_03; break;
      case  4: ud = udalt_04; break;
      case  5: ud = udalt_05; break;
      case  6: ud = udalt_06; break;
      case  7: ud = udalt_07; break;
      case  8: ud = udalt_08; break;
      case  9: ud = udalt_09; break;
      case 10: ud = udalt_10; break;
      case 11: ud = udalt_11; break;
      case 12: ud = udalt_12; break;
      case 13: ud = udalt_13; break;
      case 14: ud = udalt_14; break;
      case 15: ud = udalt_15; break;
      default: ud =  nullptr; break;
   }
   nualt = nualt + 1;
   // clang-format on
   if (nualt > 2 * maxualt)
      nualt = nualt - maxualt;

   #pragma acc parallel loop independent async\
               deviceptr(uind,ud)
   for (int i = 0; i < n; ++i) {
      ud[i][0] = uind[i][0];
      ud[i][1] = uind[i][1];
      ud[i][2] = uind[i][2];
   }
}

void ulspredSum2_acc(real (*restrict uind)[3])
{
   if (nualt < maxualt)
      return;

   constexpr double aspc[16] = {62. / 17., //
      -310. / 51.,                         //
      2170. / 323.,                        //
      -2329. / 400.,                       //
      1701. / 409.,                        //
      -806. / 323.,                        //
      1024. / 809.,                        //
      -479. / 883.,                        //
      257. / 1316.,                        //
      -434. / 7429.,                       //
      191. / 13375.,                       //
      -62. / 22287.,                       //
      3. / 7217.,                          //
      -3. / 67015.,                        //
      2. / 646323.,                        //
      -1. / 9694845.};
   constexpr double gear[6] = {6., //
      -15.,                        //
      20.,                         //
      -15.,                        //
      6.,                          //
      -1.};

   if (polpred == UPred::ASPC) {
      double c00, c01, c02, c03, c04, c05, c06, c07;
      double c08, c09, c10, c11, c12, c13, c14, c15;
      c00 = aspc[(nualt - 1 + 16) % 16];
      c01 = aspc[(nualt - 2 + 16) % 16];
      c02 = aspc[(nualt - 3 + 16) % 16];
      c03 = aspc[(nualt - 4 + 16) % 16];
      c04 = aspc[(nualt - 5 + 16) % 16];
      c05 = aspc[(nualt - 6 + 16) % 16];
      c06 = aspc[(nualt - 7 + 16) % 16];
      c07 = aspc[(nualt - 8 + 16) % 16];
      c08 = aspc[(nualt - 9 + 16) % 16];
      c09 = aspc[(nualt - 10 + 16) % 16];
      c10 = aspc[(nualt - 11 + 16) % 16];
      c11 = aspc[(nualt - 12 + 16) % 16];
      c12 = aspc[(nualt - 13 + 16) % 16];
      c13 = aspc[(nualt - 14 + 16) % 16];
      c14 = aspc[(nualt - 15 + 16) % 16];
      c15 = aspc[(nualt - 16 + 16) % 16];
      #pragma acc parallel loop independent async\
              deviceptr(uind,\
              udalt_00,udalt_01,udalt_02,udalt_03,udalt_04,udalt_05,udalt_06,\
              udalt_07,udalt_08,udalt_09,udalt_10,udalt_11,udalt_12,udalt_13,\
              udalt_14,udalt_15)
      for (int i = 0; i < n; ++i) {
         uind[i][0] = c00 * udalt_00[i][0] + c01 * udalt_01[i][0] + c02 * udalt_02[i][0] +
            c03 * udalt_03[i][0] + c04 * udalt_04[i][0] + c05 * udalt_05[i][0] +
            c06 * udalt_06[i][0] + c07 * udalt_07[i][0] + c08 * udalt_08[i][0] +
            c09 * udalt_09[i][0] + c10 * udalt_10[i][0] + c11 * udalt_11[i][0] +
            c12 * udalt_12[i][0] + c13 * udalt_13[i][0] + c14 * udalt_14[i][0] +
            c15 * udalt_15[i][0];
         uind[i][1] = c00 * udalt_00[i][1] + c01 * udalt_01[i][1] + c02 * udalt_02[i][1] +
            c03 * udalt_03[i][1] + c04 * udalt_04[i][1] + c05 * udalt_05[i][1] +
            c06 * udalt_06[i][1] + c07 * udalt_07[i][1] + c08 * udalt_08[i][1] +
            c09 * udalt_09[i][1] + c10 * udalt_10[i][1] + c11 * udalt_11[i][1] +
            c12 * udalt_12[i][1] + c13 * udalt_13[i][1] + c14 * udalt_14[i][1] +
            c15 * udalt_15[i][1];
         uind[i][2] = c00 * udalt_00[i][2] + c01 * udalt_01[i][2] + c02 * udalt_02[i][2] +
            c03 * udalt_03[i][2] + c04 * udalt_04[i][2] + c05 * udalt_05[i][2] +
            c06 * udalt_06[i][2] + c07 * udalt_07[i][2] + c08 * udalt_08[i][2] +
            c09 * udalt_09[i][2] + c10 * udalt_10[i][2] + c11 * udalt_11[i][2] +
            c12 * udalt_12[i][2] + c13 * udalt_13[i][2] + c14 * udalt_14[i][2] +
            c15 * udalt_15[i][2];
      }
   } else if (polpred == UPred::GEAR) {
      double c00, c01, c02, c03, c04, c05;
      c00 = gear[(nualt - 1 + 6) % 6];
      c01 = gear[(nualt - 2 + 6) % 6];
      c02 = gear[(nualt - 3 + 6) % 6];
      c03 = gear[(nualt - 4 + 6) % 6];
      c04 = gear[(nualt - 5 + 6) % 6];
      c05 = gear[(nualt - 6 + 6) % 6];
      #pragma acc parallel loop independent async\
              deviceptr(uind,\
              udalt_00,udalt_01,udalt_02,udalt_03,udalt_04,udalt_05)
      for (int i = 0; i < n; ++i) {
         uind[i][0] = c00 * udalt_00[i][0] + c01 * udalt_01[i][0] + c02 * udalt_02[i][0] +
            c03 * udalt_03[i][0] + c04 * udalt_04[i][0] + c05 * udalt_05[i][0];
         uind[i][1] = c00 * udalt_00[i][1] + c01 * udalt_01[i][1] + c02 * udalt_02[i][1] +
            c03 * udalt_03[i][1] + c04 * udalt_04[i][1] + c05 * udalt_05[i][1];
         uind[i][2] = c00 * udalt_00[i][2] + c01 * udalt_01[i][2] + c02 * udalt_02[i][2] +
            c03 * udalt_03[i][2] + c04 * udalt_04[i][2] + c05 * udalt_05[i][2];
      }
   } else if (polpred == UPred::LSQR) {
      using real3_ptr = real(*)[3];
      real3_ptr ppd[7] = {udalt_00, udalt_01, udalt_02, udalt_03, udalt_04, udalt_05, udalt_06};
      real3_ptr pd[7] = {ppd[(nualt - 1 + 7) % 7], ppd[(nualt - 2 + 7) % 7],
         ppd[(nualt - 3 + 7) % 7], ppd[(nualt - 4 + 7) % 7], ppd[(nualt - 5 + 7) % 7],
         ppd[(nualt - 6 + 7) % 7], ppd[(nualt - 7 + 7) % 7]};

      // k = 1 ~ 7, m = k ~ 7
      // c(k,m) = u(k) dot u(m)
      // b(1) ~ b(6) = c(1,2) ~ c(1,7)
      for (int k = 0; k < 6; ++k) {
         darray::dot(g::q0, n, &udalt_lsqr_b[k], pd[0], pd[k + 1]);
      }
      // a(1) ~ a(21)
      // OpenACC and CPU save the upper triangle.
      // c22,c23,c24,c25,c26,c27
      //     c33,c34,c35,c36,c37
      //         c44,c45,c46,c47
      //             c55,c56,c57
      //                 c66,c67
      //                     c77
      int ia = 0;
      for (int k = 1; k < 7; ++k) {
         for (int m = k; m < 7; ++m) {
            darray::dot(g::q0, n, &udalt_lsqr_a[ia], pd[k], pd[m]);
            ++ia;
         }
      }

      symlusolve<6, real>(udalt_lsqr_a, udalt_lsqr_b);

      real3_ptr pd0 = pd[0];
      real3_ptr pd1 = pd[1];
      real3_ptr pd2 = pd[2];
      real3_ptr pd3 = pd[3];
      real3_ptr pd4 = pd[4];
      real3_ptr pd5 = pd[5];
      real* bd = udalt_lsqr_b;
      #pragma acc parallel loop independent async\
                  deviceptr(uind,bd,\
                  pd0,pd1,pd2,pd3,pd4,pd5)
      for (int i = 0; i < n; ++i) {
         uind[i][0] = bd[0] * pd0[i][0] + bd[1] * pd1[i][0] + bd[2] * pd2[i][0] +
            bd[3] * pd3[i][0] + bd[4] * pd4[i][0] + bd[5] * pd5[i][0];
         uind[i][1] = bd[0] * pd0[i][1] + bd[1] * pd1[i][1] + bd[2] * pd2[i][1] +
            bd[3] * pd3[i][1] + bd[4] * pd4[i][1] + bd[5] * pd5[i][1];
         uind[i][2] = bd[0] * pd0[i][2] + bd[1] * pd1[i][2] + bd[2] * pd2[i][2] +
            bd[3] * pd3[i][2] + bd[4] * pd4[i][2] + bd[5] * pd5[i][2];
      }
   }
}
}
