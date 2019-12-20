#include "e_mplar.h"
#include "md.h"
#include "potent.h"
#include <map>
#include <tinker/detail/couple.hh>
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/polgrp.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/sizes.hh>


TINKER_NAMESPACE_BEGIN
void emplar_data(rc_op op)
{
   if (!use_potent(mpole_term) || !use_potent(polar_term))
      return;


   if (op & rc_dealloc) {
      nmdpuexclude = 0;
      device_array::deallocate(mdpuexclude, mdpuexclude_scale);
   }


   if (op & rc_alloc) {
      struct mdpu
      {
         real m, d, p, u;
      };
      auto insert_mdpu = [](std::map<std::pair<int, int>, mdpu>& a, int i,
                            int k, real val, char ch) {
         std::pair<int, int> key;
         key.first = i;
         key.second = k;
         auto it = a.find(key);
         if (it == a.end()) {
            mdpu x;
            x.m = 0;
            x.d = 0;
            x.p = 0;
            x.u = 0;
            if (ch == 'm')
               x.m = val;
            else if (ch == 'd')
               x.d = val;
            else if (ch == 'p')
               x.p = val;
            else if (ch == 'u')
               x.u = val;
            a[key] = x;
         } else {
            if (ch == 'm')
               it->second.m = val;
            else if (ch == 'd')
               it->second.d = val;
            else if (ch == 'p')
               it->second.p = val;
            else if (ch == 'u')
               it->second.u = val;
         }
      };


      // see also attach.h
      const int maxn12 = sizes::maxval;
      const int maxn13 = 3 * sizes::maxval;
      const int maxn14 = 9 * sizes::maxval;
      const int maxn15 = 27 * sizes::maxval;
      const int* couple_i12 = &couple::i12[0][0];
      const int* couple_i13 = couple::i13;
      const int* couple_i14 = couple::i14;
      const int* couple_i15 = couple::i15;
      const real m2scale = mplpot::m2scale;
      const real m3scale = mplpot::m3scale;
      const real m4scale = mplpot::m4scale;
      const real m5scale = mplpot::m5scale;


      const int maxp11 = polgrp::maxp11;
      const int maxp12 = polgrp::maxp12;
      const int maxp13 = polgrp::maxp13;
      const int maxp14 = polgrp::maxp14;
      const real d1scale = polpot::d1scale;
      const real d2scale = polpot::d2scale;
      const real d3scale = polpot::d3scale;
      const real d4scale = polpot::d4scale;
      const real p2scale = polpot::p2scale;
      const real p3scale = polpot::p3scale;
      const real p4scale = polpot::p4scale;
      const real p5scale = polpot::p5scale;
      const real p2iscale = polpot::p2iscale;
      const real p3iscale = polpot::p3iscale;
      const real p4iscale = polpot::p4iscale;
      const real p5iscale = polpot::p5iscale;
      const real u1scale = polpot::u1scale;
      const real u2scale = polpot::u2scale;
      const real u3scale = polpot::u3scale;
      const real u4scale = polpot::u4scale;


      std::map<std::pair<int, int>, mdpu> ik_scale;
      for (int i = 0; i < n; ++i) {
         int nn, bask;


         // m


         if (m2scale != 1) {
            nn = couple::n12[i];
            bask = i * maxn12;
            for (int j = 0; j < nn; ++j) {
               int k = couple_i12[bask + j];
               k -= 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, m2scale - 1, 'm');
               }
            }
         }


         if (m3scale != 1) {
            nn = couple::n13[i];
            bask = i * maxn13;
            for (int j = 0; j < nn; ++j) {
               int k = couple_i13[bask + j];
               k -= 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, m3scale - 1, 'm');
               }
            }
         }


         if (m4scale != 1) {
            nn = couple::n14[i];
            bask = i * maxn14;
            for (int j = 0; j < nn; ++j) {
               int k = couple_i14[bask + j];
               k -= 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, m4scale - 1, 'm');
               }
            }
         }


         if (m5scale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple_i15[bask + j];
               k -= 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, m5scale - 1, 'm');
               }
            }
         }


         // p


         if (p2scale != 1 || p2iscale != 1) {
            nn = couple::n12[i];
            bask = i * maxn12;
            for (int j = 0; j < nn; ++j) {
               int k = couple_i12[bask + j];
               real val = p2scale - 1;
               for (int jj = 0; jj < polgrp::np11[i]; ++jj) {
                  if (k == polgrp::ip11[i * maxp11 + jj])
                     val = p2iscale - 1;
               }
               k -= 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, val, 'p');
               }
            }
         }


         if (p3scale != 1 || p3iscale != 1) {
            nn = couple::n13[i];
            bask = i * maxn13;
            for (int j = 0; j < nn; ++j) {
               int k = couple_i13[bask + j];
               real val = p3scale - 1;
               for (int jj = 0; jj < polgrp::np11[i]; ++jj) {
                  if (k == polgrp::ip11[i * maxp11 + jj])
                     val = p3iscale - 1;
               }
               k -= 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, val, 'p');
               }
            }
         }


         if (p4scale != 1 || p4iscale != 1) {
            nn = couple::n14[i];
            bask = i * maxn14;
            for (int j = 0; j < nn; ++j) {
               int k = couple_i14[bask + j];
               real val = p4scale - 1;
               for (int jj = 0; jj < polgrp::np11[i]; ++jj) {
                  if (k == polgrp::ip11[i * maxp11 + jj])
                     val = p4iscale - 1;
               }
               k -= 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, val, 'p');
               }
            }
         }


         if (p5scale != 1 || p5iscale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple_i15[bask + j];
               real val = p5scale - 1;
               for (int jj = 0; jj < polgrp::np11[i]; ++jj) {
                  if (k == polgrp::ip11[i * maxp11 + jj])
                     val = p5iscale - 1;
               }
               k -= 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, val, 'p');
               }
            }
         }


         // d


         if (d1scale != 1) {
            nn = polgrp::np11[i];
            bask = i * maxp11;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip11[bask + j] - 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, d1scale - 1, 'd');
               }
            }
         }


         if (d2scale != 1) {
            nn = polgrp::np12[i];
            bask = i * maxp12;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip12[bask + j] - 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, d2scale - 1, 'd');
               }
            }
         }


         if (d3scale != 1) {
            nn = polgrp::np13[i];
            bask = i * maxp13;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip13[bask + j] - 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, d3scale - 1, 'd');
               }
            }
         }


         if (d4scale != 1) {
            nn = polgrp::np14[i];
            bask = i * maxp14;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip14[bask + j] - 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, d4scale - 1, 'd');
               }
            }
         }


         // u


         if (u1scale != 1) {
            nn = polgrp::np11[i];
            bask = i * maxp11;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip11[bask + j] - 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, u1scale - 1, 'u');
               }
            }
         }


         if (u2scale != 1) {
            nn = polgrp::np12[i];
            bask = i * maxp12;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip12[bask + j] - 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, u2scale - 1, 'u');
               }
            }
         }


         if (u3scale != 1) {
            nn = polgrp::np13[i];
            bask = i * maxp13;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip13[bask + j] - 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, u3scale - 1, 'u');
               }
            }
         }


         if (u4scale != 1) {
            nn = polgrp::np14[i];
            bask = i * maxp14;
            for (int j = 0; j < nn; ++j) {
               int k = polgrp::ip14[bask + j] - 1;
               if (k > i) {
                  insert_mdpu(ik_scale, i, k, u4scale - 1, 'u');
               }
            }
         }
      }
      std::vector<int> ik_vec;
      std::vector<real> scal_vec;
      for (auto& it : ik_scale) {
         ik_vec.push_back(it.first.first);
         ik_vec.push_back(it.first.second);
         scal_vec.push_back(it.second.m);
         scal_vec.push_back(it.second.d);
         scal_vec.push_back(it.second.p);
         scal_vec.push_back(it.second.u);
      }
      nmdpuexclude = ik_scale.size();
      device_array::allocate(nmdpuexclude, &mdpuexclude, &mdpuexclude_scale);
      device_array::copyin(nmdpuexclude, mdpuexclude, ik_vec.data());
      device_array::copyin(nmdpuexclude, mdpuexclude_scale, scal_vec.data());
   }
}
TINKER_NAMESPACE_END
