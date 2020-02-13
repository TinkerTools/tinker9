#include "e_polar.h"
#include "io_print.h"
#include "md.h"
#include "nblist.h"
#include "pme.h"
#include "potent.h"
#include <map>
#include <tinker/detail/couple.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/polgrp.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/sizes.hh>
#include <tinker/detail/units.hh>


TINKER_NAMESPACE_BEGIN
void epolar_data(rc_op op)
{
   if (!use_potent(polar_term))
      return;

   if (op & rc_dealloc) {
      nuexclude_ = 0;
      device_array::deallocate(uexclude_, uexclude_scale_);
      ndpexclude_ = 0;
      device_array::deallocate(dpexclude_, dpexclude_scale_);
      ndpuexclude_ = 0;
      device_array::deallocate(dpuexclude_, dpuexclude_scale_);

      device_array::deallocate(polarity, thole, pdamp, polarity_inv);

      buffer_deallocate(nep, ep, vir_ep);

      device_array::deallocate(ufld, dufld);
      device_array::deallocate(work01_, work02_, work03_, work04_, work05_,
                               work06_, work07_, work08_, work09_, work10_);
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

      struct dpu_scale
      {
         real d, p, u;
      };
      auto insert_dpu = [](std::map<std::pair<int, int>, dpu_scale>& m, int i,
                           int k, real val, char ch) {
         std::pair<int, int> key;
         key.first = i;
         key.second = k;
         auto it = m.find(key);
         if (it == m.end()) {
            dpu_scale dpu;
            dpu.d = 0;
            dpu.p = 0;
            dpu.u = 0;
            if (ch == 'd')
               dpu.d = val;
            else if (ch == 'p')
               dpu.p = val;
            else if (ch == 'u')
               dpu.u = val;
            m[key] = dpu;
         } else {
            if (ch == 'd')
               it->second.d = val;
            else if (ch == 'p')
               it->second.p = val;
            else if (ch == 'u')
               it->second.u = val;
         }
      };
      std::map<std::pair<int, int>, dpu_scale> ik_dpu;

      std::vector<int> exclik;
      std::vector<real> excls;

      u1scale = polpot::u1scale;
      u2scale = polpot::u2scale;
      u3scale = polpot::u3scale;
      u4scale = polpot::u4scale;
      exclik.clear();
      excls.clear();
      for (int i = 0; i < n; ++i) {
         int nn, bask;

         if (u1scale != 1) {
            nn = polgrp::np11[i];
            bask = i * maxp11;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip11[bask + j] - 1;
               if (k > i) {
                  insert_dpu(ik_dpu, i, k, u1scale - 1, 'u');
                  exclik.push_back(i);
                  exclik.push_back(k);
                  exclik.push_back(u1scale - 1);
               }
            }
         }

         if (u2scale != 1) {
            nn = polgrp::np12[i];
            bask = i * maxp12;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip12[bask + j] - 1;
               if (k > i) {
                  insert_dpu(ik_dpu, i, k, u2scale - 1, 'u');
                  exclik.push_back(i);
                  exclik.push_back(k);
                  exclik.push_back(u2scale - 1);
               }
            }
         }

         if (u3scale != 1) {
            nn = polgrp::np13[i];
            bask = i * maxp13;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip13[bask + j] - 1;
               if (k > i) {
                  insert_dpu(ik_dpu, i, k, u3scale - 1, 'u');
                  exclik.push_back(i);
                  exclik.push_back(k);
                  exclik.push_back(u3scale - 1);
               }
            }
         }

         if (u4scale != 1) {
            nn = polgrp::np14[i];
            bask = i * maxp14;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip14[bask + j] - 1;
               if (k > i) {
                  insert_dpu(ik_dpu, i, k, u4scale - 1, 'u');
                  exclik.push_back(i);
                  exclik.push_back(k);
                  exclik.push_back(u4scale - 1);
               }
            }
         }
      }
      nuexclude_ = excls.size();
      device_array::allocate(nuexclude_, &uexclude_, &uexclude_scale_);
      device_array::copyin(WAIT_NEW_Q, nuexclude_, uexclude_, exclik.data());
      device_array::copyin(WAIT_NEW_Q, nuexclude_, uexclude_scale_,
                           excls.data());

      d1scale = polpot::d1scale;
      d2scale = polpot::d2scale;
      d3scale = polpot::d3scale;
      d4scale = polpot::d4scale;

      p2scale = polpot::p2scale;
      p3scale = polpot::p3scale;
      p4scale = polpot::p4scale;
      p5scale = polpot::p5scale;

      p2iscale = polpot::p2iscale;
      p3iscale = polpot::p3iscale;
      p4iscale = polpot::p4iscale;
      p5iscale = polpot::p5iscale;
      exclik.clear();
      excls.clear();
      struct dp_scale
      {
         real d, p;
      };
      auto insert_dp = [](std::map<int, dp_scale>& m, int k, real val,
                          char dpchar) {
         auto it = m.find(k);
         if (it == m.end()) {
            dp_scale dp;
            dp.d = 0;
            dp.p = 0;
            if (dpchar == 'd')
               dp.d = val;
            else if (dpchar == 'p')
               dp.p = val;
            m[k] = dp;
         } else {
            if (dpchar == 'd')
               it->second.d = val;
            else if (dpchar == 'p')
               it->second.p = val;
         }
      };
      for (int i = 0; i < n; ++i) {
         std::map<int, dp_scale> k_dpscale;
         int nn, bask;

         if (d1scale != 1) {
            nn = polgrp::np11[i];
            bask = i * maxp11;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip11[bask + j] - 1;
               if (k > i) {
                  insert_dpu(ik_dpu, i, k, d1scale - 1, 'd');
                  insert_dp(k_dpscale, k, d1scale - 1, 'd');
               }
            }
         }

         if (d2scale != 1) {
            nn = polgrp::np12[i];
            bask = i * maxp12;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip12[bask + j] - 1;
               if (k > i) {
                  insert_dpu(ik_dpu, i, k, d2scale - 1, 'd');
                  insert_dp(k_dpscale, k, d2scale - 1, 'd');
               }
            }
         }

         if (d3scale != 1) {
            nn = polgrp::np13[i];
            bask = i * maxp13;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip13[bask + j] - 1;
               if (k > i) {
                  insert_dpu(ik_dpu, i, k, d3scale - 1, 'd');
                  insert_dp(k_dpscale, k, d3scale - 1, 'd');
               }
            }
         }

         if (d4scale != 1) {
            nn = polgrp::np14[i];
            bask = i * maxp14;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip14[bask + j] - 1;
               if (k > i) {
                  insert_dpu(ik_dpu, i, k, d4scale - 1, 'd');
                  insert_dp(k_dpscale, k, d4scale - 1, 'd');
               }
            }
         }

         if (p2scale != 1 || p2iscale != 1) {
            nn = couple::n12[i];
            for (int j = 0; j < nn; ++j) {
               int k = couple::i12[i][j];
               real val = p2scale - 1;
               for (int jj = 0; jj < polgrp::np11[i]; ++jj) {
                  if (k == polgrp::ip11[i * maxp11 + jj])
                     val = p2iscale - 1;
               }
               k -= 1;
               if (k > i) {
                  insert_dpu(ik_dpu, i, k, val, 'p');
                  insert_dp(k_dpscale, k, val, 'p');
               }
            }
         }

         if (p3scale != 1 || p3iscale != 1) {
            nn = couple::n13[i];
            bask = i * maxn13;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i13[bask + j];
               real val = p3scale - 1;
               for (int jj = 0; jj < polgrp::np11[i]; ++jj) {
                  if (k == polgrp::ip11[i * maxp11 + jj])
                     val = p3iscale - 1;
               }
               k -= 1;
               if (k > i) {
                  insert_dpu(ik_dpu, i, k, val, 'p');
                  insert_dp(k_dpscale, k, val, 'p');
               }
            }
         }

         if (p4scale != 1 || p4iscale != 1) {
            nn = couple::n14[i];
            bask = i * maxn14;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i14[bask + j];
               real val = p4scale - 1;
               for (int jj = 0; jj < polgrp::np11[i]; ++jj) {
                  if (k == polgrp::ip11[i * maxp11 + jj])
                     val = p4iscale - 1;
               }
               k -= 1;
               if (k > i) {
                  insert_dpu(ik_dpu, i, k, val, 'p');
                  insert_dp(k_dpscale, k, val, 'p');
               }
            }
         }

         if (p5scale != 1 || p5iscale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i15[bask + j];
               real val = p5scale - 1;
               for (int jj = 0; jj < polgrp::np11[i]; ++jj) {
                  if (k == polgrp::ip11[i * maxp11 + jj])
                     val = p5iscale - 1;
               }
               k -= 1;
               if (k > i) {
                  insert_dpu(ik_dpu, i, k, val, 'p');
                  insert_dp(k_dpscale, k, val, 'p');
               }
            }
         }

         for (auto& it : k_dpscale) {
            exclik.push_back(i);
            exclik.push_back(it.first);
            excls.push_back(it.second.d);
            excls.push_back(it.second.p);
         }
      }
      std::vector<int> dpu_ik_vec;
      std::vector<real> dpu_sc_vec;
      for (auto& it : ik_dpu) {
         dpu_ik_vec.push_back(it.first.first);
         dpu_ik_vec.push_back(it.first.second);
         dpu_sc_vec.push_back(it.second.d);
         dpu_sc_vec.push_back(it.second.p);
         dpu_sc_vec.push_back(it.second.u);
      }
      ndpuexclude_ = ik_dpu.size();
      device_array::allocate(ndpuexclude_, &dpuexclude_, &dpuexclude_scale_);
      device_array::copyin(WAIT_NEW_Q, ndpuexclude_, dpuexclude_,
                           dpu_ik_vec.data());
      device_array::copyin(WAIT_NEW_Q, ndpuexclude_, dpuexclude_scale_,
                           dpu_sc_vec.data());

      ndpexclude_ = excls.size() / 2;
      device_array::allocate(ndpexclude_, &dpexclude_, &dpexclude_scale_);
      device_array::copyin(WAIT_NEW_Q, ndpexclude_, dpexclude_, exclik.data());
      device_array::copyin(WAIT_NEW_Q, ndpexclude_, dpexclude_scale_,
                           excls.data());

      device_array::allocate(n, &polarity, &thole, &pdamp, &polarity_inv);

      buffer_allocate(&nep, &ep, &vir_ep);

      if (rc_flag & calc::grad) {
         device_array::allocate(n, &ufld, &dufld);
      } else {
         ufld = nullptr;
         dufld = nullptr;
      }

      device_array::allocate(n, &work01_, &work02_, &work03_, &work04_,
                             &work05_, &work06_, &work07_, &work08_, &work09_,
                             &work10_);
   }

   if (op & rc_init) {
      if (use_ewald()) {
         epolar_electyp = elec_t::ewald;
      } else {
         epolar_electyp = elec_t::coulomb;
      }

      udiag = polpot::udiag;

      // see also polmin in induce.f
      const double polmin = 0.00000001;
      std::vector<double> pinvbuf(n);
      for (int i = 0; i < n; ++i) {
         pinvbuf[i] = 1.0 / std::max(polar::polarity[i], polmin);
      }
      device_array::copyin(WAIT_NEW_Q, n, polarity, polar::polarity);
      device_array::copyin(WAIT_NEW_Q, n, thole, polar::thole);
      device_array::copyin(WAIT_NEW_Q, n, pdamp, polar::pdamp);
      device_array::copyin(WAIT_NEW_Q, n, polarity_inv, pinvbuf.data());
   }
}


void induce(real (*ud)[3], real (*up)[3])
{
   induce_mutual_pcg1(ud, up);

   if (inform::debug && use_potent(polar_term)) {
      std::vector<double> uindbuf;
      uindbuf.resize(3 * n);
      device_array::copyout(WAIT_NEW_Q, n, uindbuf.data(), ud);
      bool header = true;
      for (int i = 0; i < n; ++i) {
         if (polar::polarity[i] != 0) {
            if (header) {
               header = false;
               print(stdout, "\n Induced Dipole Moments (Debye) :\n");
               print(stdout,
                     "\n{0:4s}Atom{0:15s}X{0:12s}Y{0:12s}Z{0:11s}Total\n\n",
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
            print(stdout, "{:>8d}     {:13.4f}{:13.4f}{:13.4f} {:13.4f}\n",
                  i + 1, u1, u2, u3, unorm);
         }
      }
   }
}


void epolar_coulomb(int vers)
{
   extern void epolar_coulomb_acc(int, const real(*)[3], const real(*)[3]);


   induce(uind, uinp);
#if TINKER_CUDART
   if (mlist_version() == NBList::spatial) {
      extern void epolar_coulomb_cu(int, const real(*)[3], const real(*)[3]);
      epolar_coulomb_cu(vers, uind, uinp);
   } else
#endif
      epolar_coulomb_acc(vers, uind, uinp);
}


void epolar(int vers)
{
   if (epolar_electyp == elec_t::coulomb)
      epolar_coulomb(vers);
   else if (epolar_electyp == elec_t::ewald)
      epolar_ewald(vers);
}
TINKER_NAMESPACE_END
