#include "ff/hippo/edisp.h"
#include "ff/energy.h"
#include "ff/pme.h"
#include "ff/potent.h"
#include "math/zero.h"
#include "tool/externfunc.h"
#include <cassert>
#include <tinker/detail/couple.hh>
#include <tinker/detail/disp.hh>
#include <tinker/detail/dsppot.hh>
#include <tinker/detail/limits.hh>
#include <tinker/detail/mutant.hh>
#include <tinker/routines.h>

namespace tinker {
// The C6 sums of the environment and of the mutated ligand, and the copy of the
// forward transform that separates their two transforms (edisp1.f:edrecip1).
static double dsp_sume, dsp_suml;
static real* dsp_qfwd;

bool useDEwald()
{
   return use(Potent::DISP) and static_cast<bool>(limits::use_dewald);
}

void edispData(RcOp op)
{
   if (not use(Potent::DISP))
      return;

   auto rc_a = rc_flag & calc::analyz;

   if (op & RcOp::DEALLOC) {
      darray::deallocate(csix, adisp);
      darray::deallocate(dspexclude, dspexclude_scale);
      darray::deallocate(dsp_qfwd);
      dsp_qfwd = nullptr;

      if (rc_a) {
         bufferDeallocate(rc_flag, ndisp);
         bufferDeallocate(rc_flag, edsp, vir_edsp, dedspx, dedspy, dedspz);
      }
      ndisp = nullptr;
      edsp = nullptr;
      vir_edsp = nullptr;
      dedspx = nullptr;
      dedspy = nullptr;
      dedspz = nullptr;

      elrc_vol_dsp = 0;
      vlrc_vol_dsp = 0;
   }

   if (op & RcOp::ALLOC) {
      darray::allocate(n, &csix, &adisp);
      // a mutated ligand may be weighted apart from the environment (dispSplit)
      dsp_qfwd = nullptr;
      if (useDEwald() and mutant::nmut > 0) {
         const auto& st = *dpme_unit;
         darray::allocate(2 * st.nfft1 * st.nfft2 * st.nfft3, &dsp_qfwd);
      }

      ndisp = nullptr;
      edsp = eng_buf_vdw;
      vir_edsp = vir_buf_vdw;
      dedspx = gx_vdw;
      dedspy = gy_vdw;
      dedspz = gz_vdw;
      if (rc_a) {
         bufferAllocate(rc_flag, &ndisp);
         bufferAllocate(rc_flag, &edsp, &vir_edsp, &dedspx, &dedspy, &dedspz);
      }

      if (dsppot::use_dcorr && not useDEwald()) {
         double elrc = 0, vlrc = 0, delrc = 0, dvlrc = 0;
         tinker_f_evcorr1({const_cast<char*>("DISP"), 4}, &elrc, &vlrc, &delrc, &dvlrc);
         elrc_vol_dsp = elrc * boxVolume();
         vlrc_vol_dsp = vlrc * boxVolume();
      } else {
         elrc_vol_dsp = 0;
         vlrc_vol_dsp = 0;
      }
      dsp2scale = dsppot::dsp2scale;
      dsp3scale = dsppot::dsp3scale;
      dsp4scale = dsppot::dsp4scale;
      dsp5scale = dsppot::dsp5scale;
      std::vector<int> exclik;
      std::vector<real> excls;
      // see also attach.f
      const int maxn13 = 3 * sizes::maxval;
      const int maxn14 = 9 * sizes::maxval;
      const int maxn15 = 27 * sizes::maxval;
      for (int i = 0; i < n; ++i) {
         int nn;
         int bask;

         if (dsp2scale != 1) {
            nn = couple::n12[i];
            for (int j = 0; j < nn; ++j) {
               int k = couple::i12[i][j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(dsp2scale);
               }
            }
         }

         if (dsp3scale != 1) {
            nn = couple::n13[i];
            bask = i * maxn13;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i13[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(dsp3scale);
               }
            }
         }

         if (dsp4scale != 1) {
            nn = couple::n14[i];
            bask = i * maxn14;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i14[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(dsp4scale);
               }
            }
         }

         if (dsp5scale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i15[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(dsp5scale);
               }
            }
         }
      }
      ndspexclude = excls.size();
      darray::allocate(ndspexclude, &dspexclude, &dspexclude_scale);
      darray::copyin(g::q0, ndspexclude, dspexclude, exclik.data());
      darray::copyin(g::q0, ndspexclude, dspexclude_scale, excls.data());
      waitFor(g::q0);
   }

   if (op & RcOp::INIT) {
      csixpr = disp::csixpr;
      darray::copyin(g::q0, n, csix, disp::csix);
      darray::copyin(g::q0, n, adisp, disp::adisp);
      waitFor(g::q0);

      dsp_sume = 0;
      dsp_suml = 0;
      if (mutant::nmut > 0) {
         for (int i = 0; i < n; ++i) {
            if (mutant::mutg[i] == 0)
               dsp_sume += disp::csix[i];
            else
               dsp_suml += disp::csix[i];
         }
      }
   }
}
}

namespace tinker {
TINKER_FVOID2(acc1, cu1, edispNonEwald, int);
static void edispNonEwald(int vers)
{
   TINKER_FCALL2(acc1, cu1, edispNonEwald, vers);
}

// The reciprocal and self terms weight the ligand by the soft-core weight v of
// the real space pairs: v between the ligand and the environment, and within
// the ligand v when it is annihilated or 1 when it is decoupled. The weighted
// C6 sum product replaces csixpr (edisp.f:dspmutwt).
static DispSplit dispSplit(real& cpr)
{
   DispSplit ds = {false, 1, 1, nullptr, nullptr};
   cpr = csixpr;
   if (mutant::nmut == 0)
      return ds;
   double vl = vlam;
   double v = vl * vl * vl * vl / std::sqrt(1 + vl * vl - vl * vl * vl);
   if (v == 1)
      return ds;

   double wel = v;
   double wll = (vcouple == Vdw::ANNIHILATE ? v : 1);
   assert(dsp_qfwd);
   ds.split = true;
   ds.wel = wel;
   ds.wll = wll;
   ds.part = mut;
   ds.qfwd = reinterpret_cast<const real(*)[2]>(dsp_qfwd);
   cpr = dsp_sume * dsp_sume + 2 * wel * dsp_sume * dsp_suml + wll * dsp_suml * dsp_suml;
   return ds;
}

TINKER_FVOID2(acc1, cu1, edispEwaldReal, int);
TINKER_FVOID2(acc1, cu1, edispEwaldRecipSelf, int, DispSplit);
TINKER_FVOID2(acc1, cu1, pmeConvDisp, int, DispSplit);
static void edispEwald(int vers)
{
   TINKER_FCALL2(acc1, cu1, edispEwaldReal, vers);

   // recip and self
   auto do_e = vers & calc::energy;
   auto do_v = vers & calc::virial;
   auto do_g = vers & calc::grad;
   PMEUnit u = dpme_unit;

   real cpr;
   const DispSplit ds = dispSplit(cpr);
   gridDisp(u, csix, ds.part);
   fftfront(u);
   // each point of the convolution reads the opposite wavevector as well
   if (ds.split)
      darray::copy(g::q0, 2 * u->nfft1 * u->nfft2 * u->nfft3, dsp_qfwd, u->qgrid);
   TINKER_FCALL2(acc1, cu1, pmeConvDisp, vers, ds);
   if (do_g) {
      fftback(u);
   }
   TINKER_FCALL2(acc1, cu1, edispEwaldRecipSelf, vers, ds);

   // account for the total energy and virial correction term
   if CONSTEXPR (do_e || do_v) {
      const real aewald = u->aewald;
      const real denom0 = 6 * boxVolume() / std::pow(M_PI, 1.5);
      energy_prec term = cpr * aewald * aewald * aewald / denom0;
      if CONSTEXPR (do_e) {
         energy_edsp -= term;
         energy_vdw -= term;
      }
      if CONSTEXPR (do_v) {
         virial_edsp[0] += term; // xx
         virial_edsp[4] += term; // yy
         virial_edsp[8] += term; // zz
         virial_vdw[0] += term;  // xx
         virial_vdw[4] += term;  // yy
         virial_vdw[8] += term;  // zz
      }
   }
}

void edisp(int vers)
{
   auto rc_a = rc_flag & calc::analyz;
   auto do_a = vers & calc::analyz;
   auto do_e = vers & calc::energy;
   auto do_v = vers & calc::virial;
   auto do_g = vers & calc::grad;

   zeroOnHost(energy_edsp, virial_edsp);
   size_t bsize = bufferSize();
   if (rc_a) {
      if (do_a)
         darray::zero(g::q0, bsize, ndisp);
      if (do_e)
         darray::zero(g::q0, bsize, edsp);
      if (do_v)
         darray::zero(g::q0, bsize, vir_edsp);
      if (do_g)
         darray::zero(g::q0, n, dedspx, dedspy, dedspz);
   }

   if (useDEwald())
      edispEwald(vers);
   else
      edispNonEwald(vers);

   if (do_e) {
      if (elrc_vol_dsp != 0) {
         energy_prec corr = elrc_vol_dsp / boxVolume();
         energy_edsp += corr;
         energy_vdw += corr;
      }
   }
   if (do_v) {
      if (vlrc_vol_dsp != 0) {
         virial_prec term = vlrc_vol_dsp / boxVolume();
         virial_edsp[0] += term; // xx
         virial_edsp[4] += term; // yy
         virial_edsp[8] += term; // zz
         virial_vdw[0] += term;
         virial_vdw[4] += term;
         virial_vdw[8] += term;
      }
   }
   if (rc_a) {
      if (do_e) {
         EnergyBuffer u = edsp;
         energy_prec e = energyReduce(u);
         energy_edsp += e;
         energy_vdw += e;
      }
      if (do_v) {
         VirialBuffer u = vir_edsp;
         virial_prec v[9];
         virialReduce(v, u);
         for (int iv = 0; iv < 9; ++iv) {
            virial_edsp[iv] += v[iv];
            virial_vdw[iv] += v[iv];
         }
      }
      if (do_g)
         sumGradient(gx_vdw, gy_vdw, gz_vdw, dedspx, dedspy, dedspz);
   }
}
}
