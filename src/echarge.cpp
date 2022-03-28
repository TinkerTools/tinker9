#include "ff/pchg/echarge.h"
#include "ff/elec.h"
#include "ff/energy.h"
#include "ff/nblist.h"
#include "ff/pchg/echglj.h"
#include "ff/pme.h"
#include "ff/potent.h"
#include "tool/zero.h"
#include <tinker/detail/chgpot.hh>
#include <tinker/detail/couple.hh>
#include <tinker/detail/sizes.hh>

namespace tinker {
void echargeData(RcOp op)
{
   if (!usePotent(Potent::CHARGE))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      ncexclude = 0;
      darray::deallocate(cexclude, cexclude_scale);

      if (rc_a) {
         buffer_deallocate(rc_flag, nec);
         buffer_deallocate(rc_flag, ec, vir_ec, decx, decy, decz);
      }
      nec = nullptr;
      ec = nullptr;
      vir_ec = nullptr;
      decx = nullptr;
      decy = nullptr;
      decz = nullptr;
   }

   if (op & rc_alloc) {
      ebuffer = chgpot::ebuffer;

      c2scale = chgpot::c2scale;
      c3scale = chgpot::c3scale;
      c4scale = chgpot::c4scale;
      c5scale = chgpot::c5scale;

      std::vector<int> exclik;
      std::vector<real> excl;
      // see attach.f
      const int maxn13 = 3 * sizes::maxval;
      const int maxn14 = 9 * sizes::maxval;
      const int maxn15 = 27 * sizes::maxval;
      for (int i = 0; i < n; ++i) {
         int nn;
         int bask;
         if (c2scale != 1) {
            nn = couple::n12[i];
            for (int j = 0; j < nn; ++j) {
               int k = couple::i12[i][j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excl.push_back(c2scale);
               }
            }
         }
         if (c3scale != 1) {
            nn = couple::n13[i];
            bask = i * maxn13;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i13[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excl.push_back(c3scale);
               }
            }
         }
         if (c4scale != 1) {
            nn = couple::n14[i];
            bask = i * maxn14;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i14[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excl.push_back(c4scale);
               }
            }
         }
         if (c5scale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i15[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excl.push_back(c5scale);
               }
            }
         }
      }
      ncexclude = excl.size();
      darray::allocate(ncexclude, &cexclude, &cexclude_scale);
      darray::copyin(g::q0, ncexclude, cexclude, exclik.data());
      darray::copyin(g::q0, ncexclude, cexclude_scale, excl.data());
      wait_for(g::q0);

      nec = nullptr;
      ec = eng_buf_elec;
      vir_ec = vir_buf_elec;
      decx = gx_elec;
      decy = gy_elec;
      decz = gz_elec;
      if (rc_a) {
         buffer_allocate(rc_flag, &nec);
         buffer_allocate(rc_flag, &ec, &vir_ec, &decx, &decy, &decz);
      }
   }

   if (op & rc_init) {}
}

void echarge(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_a = vers & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   zeroOnHost(energy_ec, virial_ec);
   size_t bsize = buffer_size();
   if (rc_a) {
      if (do_a)
         darray::zero(g::q0, bsize, nec);
      if (do_e)
         darray::zero(g::q0, bsize, ec);
      if (do_v)
         darray::zero(g::q0, bsize, vir_ec);
      if (do_g)
         darray::zero(g::q0, n, decx, decy, decz);
   }

   if (useEwald()) {
      echarge_ewald_recip_self(vers);
#if TINKER_CUDART
      if (clistVersion() & Nbl::SPATIAL)
         echarge_ewald_real_cu(vers);
      else
#endif
         echarge_ewald_real_acc(vers);
      pme_stream_finish_wait(use_pme_stream and (vers & calc::analyz));
   } else
      echarge_nonewald(vers);

   if (rc_a) {
      if (do_e) {
         energy_buffer u = ec;
         energy_prec e = energy_reduce(u);
         energy_ec += e;
         energy_elec += e;
      }
      if (do_v) {
         virial_buffer u = vir_ec;
         virial_prec v[9];
         virial_reduce(v, u);
         for (int iv = 0; iv < 9; ++iv) {
            virial_ec[iv] += v[iv];
            virial_elec[iv] += v[iv];
         }
      }
      if (do_g)
         sum_gradient(gx_elec, gy_elec, gz_elec, decx, decy, decz);
   }
}

void echarge_nonewald(int vers)
{
#if TINKER_CUDART
   if (clistVersion() & Nbl::SPATIAL)
      echarge_nonewald_cu(vers);
   else
#endif
      echarge_nonewald_acc(vers);
}

void echarge_ewald_recip_self(int vers)
{
   pme_stream_start_wait(use_pme_stream);

   // ewald recip space, self term
   // ewald real space

   const PMEUnit pu = epme_unit;
   grid_pchg(pu, pchg);
   fftfront(pu);
   if (vers & calc::virial) {
      if (vers & calc::energy) {
         pme_conv(pu, ec, vir_ec);
      } else {
         pme_conv(pu, vir_ec);
      }
   } else {
      if (vers & calc::energy) {
         pme_conv(pu, ec);
      } else {
         pme_conv(pu);
      }
   }
   if (vers & calc::grad) {
      fftback(pu);
   }

   // fphi_pchg, recip, self

#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      echarge_ewald_fphi_self_cu(vers);
   else
#endif
      echarge_ewald_fphi_self_acc(vers);

   pme_stream_finish_record(use_pme_stream);
}
}
