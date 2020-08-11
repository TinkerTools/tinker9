#include "epolar.h"
#include "empole.h"
#include "md.h"
#include "nblist.h"
#include "pme.h"
#include "potent.h"
#include "tool/host_zero.h"
#include "tool/io_print.h"
#include <map>
#include <tinker/detail/couple.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/polgrp.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/sizes.hh>
#include <tinker/detail/units.hh>


namespace tinker {
void epolar_data(rc_op op)
{
   if (!use_potent(polar_term))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      ndwexclude = 0;
      darray::deallocate(dwexclude, dwexclude_scale);
      
      // darray::deallocate(polarity, thole, pdamp, polarity_inv);

      if (rc_a) {
         buffer_deallocate(rc_flag, nep);
         buffer_deallocate(rc_flag, ep, vir_ep, depx, depy, depz);
      }
      nep = nullptr;
      ep = nullptr;
      vir_ep = nullptr;
      depx = nullptr;
      depy = nullptr;
      depz = nullptr;

      darray::deallocate(ufld, dufld);
      darray::deallocate(work01_, work02_, work03_, work04_, work05_, work06_,
                         work07_, work08_, work09_, work10_);
   }

   if (op & rc_alloc) {
      // see also attach.h
      const int maxn13 = 3 * sizes::maxval;
      const int maxn14 = 9 * sizes::maxval;
      const int maxn15 = 27 * sizes::maxval;
      const int maxp11 = polgrp::maxp11;
      const int maxp12 = polgrp::maxp12;
      const int maxp13 = polgrp::maxp13;
      const int maxp14 = polgrp::maxp14;

      
      d1scale = polpot::d1scale;
      d2scale = polpot::d2scale;
      d3scale = polpot::d3scale;
      d4scale = polpot::d4scale;

      w2scale = polpot::w2scale;
      w3scale = polpot::w3scale;
      w4scale = polpot::w4scale;
      w5scale = polpot::w5scale;

      exclik.clear();
      excls.clear();
      struct dw_scale
      {
         real d, w;
      };
      auto insert_dw = [](std::map<int, dw_scale>& m, int k, real val,
                          char dwchar) {
         auto it = m.find(k);
         if (it == m.end()) {
            dw_scale dw;
            dw.d = 0;
            dw.w = 0;
            if (dwchar == 'd')
               dw.d = val;
            else if (dwchar == 'w')
               dw.w = val;
            m[k] = dw;
         } else {
            if (dwchar == 'd')
               it->second.d = val;
            else if (dwchar == 'w')
               it->second.p = val;
         }
      };
      for (int i = 0; i < n; ++i) {
         std::map<int, dw_scale> k_dwscale;
         int nn, bask;

         if (d1scale != 1) {
            nn = polgrp::np11[i];
            bask = i * maxp11;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip11[bask + j] - 1;
               if (k > i) 
                  insert_dw(k_dwscale, k, d1scale - 1, 'd');
            }
         }

         if (d2scale != 1) {
            nn = polgrp::np12[i];
            bask = i * maxp12;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip12[bask + j] - 1;
               if (k > i) 
                  insert_dw(k_dwscale, k, d2scale - 1, 'd');
            }
         }

         if (d3scale != 1) {
            nn = polgrp::np13[i];
            bask = i * maxp13;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip13[bask + j] - 1;
               if (k > i) 
                  insert_dw(k_dwscale, k, d3scale - 1, 'd');
            }
         }

         if (d4scale != 1) {
            nn = polgrp::np14[i];
            bask = i * maxp14;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip14[bask + j] - 1;
               if (k > i) 
                  insert_dw(k_dwscale, k, d4scale - 1, 'd');
            }
         }

         if (w2scale != 1) {
            nn = couple::n12[i];
            bask = i * maxp12;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i12[bask + j] - 1;
               if (k > i) 
                  insert_dw(k_dwscale, k, w2scale - 1, 'w');
            }
         }

         if (w3scale != 1) {
            nn = couple::n13[i];
            bask = i * maxp13;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i13[bask + j] - 1;
               if (k > i) 
                  insert_dw(k_dwscale, k, w3scale - 1, 'w');
            }
         }

         if (w4scale != 1) {
            nn = couple::n14[i];
            bask = i * maxp14;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i14[bask + j] - 1;
               if (k > i) 
                  insert_dw(k_dwscale, k, w4scale - 1, 'w');
            }
         }

         if (w5scale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i15[bask + j] - 1;
               if (k > i) 
                  insert_dw(k_dwscale, k, w5scale - 1, 'w');
            }
         }

         

         for (auto& it : k_dwscale) {
            exclik.push_back(i);
            exclik.push_back(it.first);
            excls.push_back(it.second.d);
            excls.push_back(it.second.w);
         }
      }
      

      ndwexclude = excls.size() / 2;
      darray::allocate(ndwexclude, &dwexclude, &dwexclude_scale);
      darray::copyin(WAIT_NEW_Q, ndwexclude, dwexclude, exclik.data());
      darray::copyin(WAIT_NEW_Q, ndwexclude, dwexclude_scale, excls.data());

      darray::allocate(n, &polarity, &thole, &pdamp, &polarity_inv);

      nep = nullptr;
      ep = eng_buf_elec;
      vir_ep = vir_buf_elec;
      depx = gx_elec;
      depy = gy_elec;
      depz = gz_elec;
      if (rc_a) {
         buffer_allocate(rc_flag, &nep);
         buffer_allocate(rc_flag, &ep, &vir_ep, &depx, &depy, &depz);
      }

      if (rc_flag & calc::grad) {
         darray::allocate(n, &ufld, &dufld);
      } else {
         ufld = nullptr;
         dufld = nullptr;
      }

      darray::allocate(n, &work01_, &work02_, &work03_, &work04_, &work05_,
                       &work06_, &work07_, &work08_, &work09_, &work10_);
   }

   if (op & rc_init) {
      udiag = polpot::udiag;

      const double polmin = 1.0e-16;
      std::vector<double> pinvbuf(n);
      for (int i = 0; i < n; ++i) {
         pinvbuf[i] = 1.0 / std::max(polar::polarity[i], polmin);
      }
      darray::copyin(WAIT_NEW_Q, n, polarity, polar::polarity);
      darray::copyin(WAIT_NEW_Q, n, thole, polar::thole);
      darray::copyin(WAIT_NEW_Q, n, pdamp, polar::pdamp);
      darray::copyin(WAIT_NEW_Q, n, polarity_inv, pinvbuf.data());
   }
}


void induce(real (*ud)[3], real (*up)[3])
{
   induce_mutual_pcg1(ud, up);

   if (inform::debug && use_potent(polar_term)) {
      std::vector<double> uindbuf;
      uindbuf.resize(3 * n);
      darray::copyout(WAIT_NEW_Q, n, uindbuf.data(), ud);
      bool header = true;
      for (int i = 0; i < n; ++i) {
         if (polar::polarity[i] != 0) {
            if (header) {
               header = false;
               print(stdout, "\n Induced Dipole Moments (Debye) :\n");
               print(stdout,
                     "\n    Atom %1$13s X %1$10s Y %1$10s Z %1$9s Total\n\n",
                     "");
            }
            double u1 = uindbuf[3 * i];
            double u2 = uindbuf[3 * i + 1];
            double u3 = uindbuf[3 * i + 2];
            double unorm = std::sqrt(u1 * u1 + u2 * u2 + u3 * u3);
            u1 *= units::debye;
            u2 *= units::debye;
            u3 *= units::debye;
            unorm *= units::debye;
            print(stdout, "%8d     %13.4f%13.4f%13.4f %13.4f\n", i + 1, u1, u2,
                  u3, unorm);
         }
      }
   }
}


void epolar(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_a = vers & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   host_zero(energy_ep, virial_ep);
   size_t bsize = buffer_size();
   if (rc_a) {
      if (do_a)
         darray::zero(PROCEED_NEW_Q, bsize, nep);
      if (do_e)
         darray::zero(PROCEED_NEW_Q, bsize, ep);
      if (do_v) {
         darray::zero(PROCEED_NEW_Q, bsize, vir_ep);
      }
      if (do_g) {
         darray::zero(PROCEED_NEW_Q, n, depx, depy, depz);
      }
   }


   mpole_init(vers);
   if (use_ewald())
      epolar_ewald(vers);
   else
      epolar_nonewald(vers);
   torque(vers, depx, depy, depz);
   if (do_v) {
      virial_buffer u2 = vir_trq;
      virial_prec v2[9];
      virial_reduce(v2, u2);
      for (int iv = 0; iv < 9; ++iv) {
         virial_ep[iv] += v2[iv];
         virial_elec[iv] += v2[iv];
      }
   }


   if (rc_a) {
      if (do_e) {
         energy_buffer u = ep;
         energy_prec e = energy_reduce(u);
         energy_ep += e;
         energy_elec += e;
      }
      if (do_v) {
         virial_buffer u1 = vir_ep;
         virial_prec v1[9];
         virial_reduce(v1, u1);
         for (int iv = 0; iv < 9; ++iv) {
            virial_ep[iv] = v1[iv];
            virial_elec[iv] += v1[iv];
         }
      }
      if (do_g)
         sum_gradient(gx_elec, gy_elec, gz_elec, depx, depy, depz);
   }
}


void epolar_nonewald(int vers)
{
   // v0: E_dot
   // v1: EGV = E_dot + GV
   // v3: EA = E_pair + A
   // v4: EG = E_dot + G
   // v5: G
   // v6: GV
   bool edot = vers & calc::energy; // if not do_e, edot = false
   if (vers & calc::energy && vers & calc::analyz)
      edot = false; // if do_e and do_a, edot = false
   int ver2 = vers;
   if (edot)
      ver2 &= ~calc::energy; // toggle off the calc::energy flag

   induce(uind, uinp);
   if (edot)
      epolar0_dotprod(uind, udirp);
   if (vers != calc::v0) {
#if TINKER_CUDART
      if (mlist_version() & NBL_SPATIAL)
         epolar_nonewald_cu(ver2, uind, uinp);
      else
#endif
         epolar_nonewald_acc(ver2, uind, uinp);
   }
}


void epolar_ewald(int vers)
{
   // v0: E_dot
   // v1: EGV = E_dot + GV
   // v3: EA = E_pair + A
   // v4: EG = E_dot + G
   // v5: G
   // v6: GV
   bool edot = vers & calc::energy; // if not do_e, edot = false
   if (vers & calc::energy && vers & calc::analyz)
      edot = false; // if do_e and do_a, edot = false
   int ver2 = vers;
   if (edot)
      ver2 &= ~calc::energy; // toggle off the calc::energy flag

   induce(uind, uinp);
   if (edot)
      epolar0_dotprod(uind, udirp);
   if (vers != calc::v0) {
      epolar_ewald_real(ver2);
      epolar_ewald_recip_self(ver2);
   }
}


void epolar_ewald_real(int vers)
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      epolar_ewald_real_cu(vers, uind, uinp);
   else
#endif
      epolar_ewald_real_acc(vers, uind, uinp);
}


void epolar_ewald_recip_self(int vers)
{
   epolar_ewald_recip_self_acc(vers, uind, uinp);
}


void epolar0_dotprod(const real (*uind)[3], const real (*udirp)[3])
{
   return epolar0_dotprod_acc(uind, udirp);
}
}
