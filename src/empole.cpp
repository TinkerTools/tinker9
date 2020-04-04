#include "empole.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"
#include <tinker/detail/couple.hh>
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/sizes.hh>


TINKER_NAMESPACE_BEGIN
void empole_data(rc_op op)
{
   if (!use_potent(mpole_term))
      return;

   if (op & rc_dealloc) {
      nmexclude = 0;
      darray::deallocate(mexclude, mexclude_scale);

      buffer_deallocate(nem, em, vir_em);
   }

   if (op & rc_alloc) {
      m2scale = mplpot::m2scale;
      m3scale = mplpot::m3scale;
      m4scale = mplpot::m4scale;
      m5scale = mplpot::m5scale;

      std::vector<int> exclik;
      std::vector<real> excls;
      // see also attach.f
      const int maxn13 = 3 * sizes::maxval;
      const int maxn14 = 9 * sizes::maxval;
      const int maxn15 = 27 * sizes::maxval;
      for (int i = 0; i < n; ++i) {
         int nn;
         int bask;

         if (m2scale != 1) {
            nn = couple::n12[i];
            for (int j = 0; j < nn; ++j) {
               int k = couple::i12[i][j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(m2scale - 1);
               }
            }
         }

         if (m3scale != 1) {
            nn = couple::n13[i];
            bask = i * maxn13;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i13[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(m3scale - 1);
               }
            }
         }

         if (m4scale != 1) {
            nn = couple::n14[i];
            bask = i * maxn14;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i14[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(m4scale - 1);
               }
            }
         }

         if (m5scale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i15[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(m5scale - 1);
               }
            }
         }
      }
      nmexclude = excls.size();
      darray::allocate(nmexclude, &mexclude, &mexclude_scale);
      darray::copyin(WAIT_NEW_Q, nmexclude, mexclude, exclik.data());
      darray::copyin(WAIT_NEW_Q, nmexclude, mexclude_scale, excls.data());

      buffer_allocate(&nem, &em, &vir_em);
   }

   if (op & rc_init) {
   }
}


void empole(int vers)
{
   if (use_ewald())
      empole_ewald(vers);
   else
      empole_nonewald(vers);
}


void empole_nonewald(int vers)
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      empole_nonewald_cu(vers);
   else
#endif
      empole_nonewald_acc(vers);
}


void empole_ewald(int vers)
{
   empole_ewald_real_self(vers);
   empole_ewald_recip(vers);
}


void empole_ewald_real_self(int vers)
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      empole_ewald_real_self_cu(vers);
   else
#endif
      empole_ewald_real_self_acc(vers);
}


void empole_ewald_recip(int vers)
{
   empole_ewald_recip_acc(vers);
}
TINKER_NAMESPACE_END
