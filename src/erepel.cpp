#include "ff/hippo/erepel.h"
#include "ff/amoeba/empole.h"
#include "ff/elec.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "md/inc.h"
#include "mod/repel.h"
#include "tool/darray.h"
#include "tool/zero.h"
#include <tinker/detail/couple.hh>
#include <tinker/detail/repel.hh>
#include <tinker/detail/reppot.hh>
#include <tinker/detail/sizes.hh>

namespace tinker {
void erepelData(RcOp op)
{
   if (!use_potent(repuls_term))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      darray::deallocate(sizpr, dmppr, elepr);
      nrepexclude = 0;
      darray::deallocate(repexclude, repexclude_scale);

      if (rc_a) {
         buffer_deallocate(rc_flag, nrep);
         buffer_deallocate(rc_flag, er, vir_er, derx, dery, derz);
      }
      nrep = nullptr;
      er = nullptr;
      vir_er = nullptr;
      derx = nullptr;
      dery = nullptr;
      derz = nullptr;
   }

   if (op & rc_alloc) {
      darray::allocate(n, &sizpr, &dmppr, &elepr);

      nrep = nullptr;
      er = eng_buf_vdw;
      vir_er = vir_buf_vdw;
      derx = gx_vdw;
      dery = gy_vdw;
      derz = gz_vdw;
      if (rc_a) {
         buffer_allocate(rc_flag, &nrep);
         buffer_allocate(rc_flag, &er, &vir_er, &derx, &dery, &derz);
      }

      r2scale = reppot::r2scale;
      r3scale = reppot::r3scale;
      r4scale = reppot::r4scale;
      r5scale = reppot::r5scale;
      std::vector<int> exclik;
      std::vector<real> excls;
      // see also attach.f
      const int maxn13 = 3 * sizes::maxval;
      const int maxn14 = 9 * sizes::maxval;
      const int maxn15 = 27 * sizes::maxval;
      for (int i = 0; i < n; ++i) {
         int nn;
         int bask;

         if (r2scale != 1) {
            nn = couple::n12[i];
            for (int j = 0; j < nn; ++j) {
               int k = couple::i12[i][j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(r2scale);
               }
            }
         }

         if (r3scale != 1) {
            nn = couple::n13[i];
            bask = i * maxn13;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i13[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(r3scale);
               }
            }
         }

         if (r4scale != 1) {
            nn = couple::n14[i];
            bask = i * maxn14;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i14[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(r4scale);
               }
            }
         }

         if (r5scale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i15[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(r5scale);
               }
            }
         }
      }
      nrepexclude = excls.size();
      darray::allocate(nrepexclude, &repexclude, &repexclude_scale);
      darray::copyin(g::q0, nrepexclude, repexclude, exclik.data());
      darray::copyin(g::q0, nrepexclude, repexclude_scale, excls.data());
      wait_for(g::q0);
   }

   if (op & rc_init) {
      darray::copyin(g::q0, n, sizpr, repel::sizpr);
      darray::copyin(g::q0, n, dmppr, repel::dmppr);
      darray::copyin(g::q0, n, elepr, repel::elepr);
      wait_for(g::q0);
   }
}

void erepel_cu(int);
void erepel_acc(int);
void erepel(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_a = vers & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   zeroOnHost(energy_er, virial_er);
   size_t bsize = buffer_size();
   if (rc_a) {
      if (do_a)
         darray::zero(g::q0, bsize, nrep);
      if (do_e)
         darray::zero(g::q0, bsize, er);
      if (do_v)
         darray::zero(g::q0, bsize, vir_er);
      if (do_g)
         darray::zero(g::q0, n, derx, dery, derz);
   }

   mpoleInit(vers);
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      erepel_cu(vers);
   else
#endif
      erepel_acc(vers);

   torque(vers, derx, dery, derz);
   if (do_v) {
      virial_buffer u2 = vir_trq;
      virial_prec v2[9];
      virial_reduce(v2, u2);
      for (int iv = 0; iv < 9; ++iv) {
         virial_er[iv] += v2[iv];
         virial_vdw[iv] += v2[iv];
      }
   }

   if (rc_a) {
      if (do_e) {
         energy_buffer u = er;
         energy_prec e = energy_reduce(u);
         energy_er += e;
         energy_vdw += e;
      }
      if (do_v) {
         virial_buffer u = vir_er;
         virial_prec v[9];
         virial_reduce(v, u);
         for (int iv = 0; iv < 9; ++iv) {
            virial_er[iv] += v[iv];
            virial_vdw[iv] += v[iv];
         }
      }
      if (do_g)
         sum_gradient(gx_vdw, gy_vdw, gz_vdw, derx, dery, derz);
   }
}
}
