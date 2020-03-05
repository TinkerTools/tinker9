#include "e_vdw.h"
#include "io_fort_str.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"
#include "tinker_rt.h"
#include <cassert>
#include <map>
#include <tinker/detail/couple.hh>
#include <tinker/detail/mutant.hh>
#include <tinker/detail/sizes.hh>
#include <tinker/detail/vdw.hh>
#include <tinker/detail/vdwpot.hh>

TINKER_NAMESPACE_BEGIN
void evdw_data(rc_op op)
{
   if (!use_potent(vdw_term))
      return;

   typedef int new_type;
   typedef int old_type;
   static std::map<old_type, new_type> jmap;
   static std::vector<new_type> jvec;
   static std::vector<new_type> jvdwbuf;
   static int jcount;

   if (op & rc_dealloc) {
      // local static members
      jmap.clear();
      jvec.clear();
      jvdwbuf.clear();
      jcount = 0;

      darray::deallocate(ired, kred, xred, yred, zred, gxred, gyred, gzred,
                         jvdw, radmin, epsilon, vlam);

      nvexclude_ = 0;
      darray::deallocate(vexclude_, vexclude_scale_);

      buffer_deallocate(nev, ev, vir_ev);

      elrc_vol = 0;
      vlrc_vol = 0;
   }

   if (op & rc_alloc) {
      darray::allocate(n, &ired, &kred, &xred, &yred, &zred);
      if (rc_flag & calc::grad) {
         darray::allocate(n, &gxred, &gyred, &gzred);
      } else {
         gxred = nullptr;
         gyred = nullptr;
         gzred = nullptr;
      }

      darray::allocate(n, &jvdw);

      jvdwbuf.resize(n);
      assert(jmap.size() == 0);
      assert(jvec.size() == 0);
      jcount = 0;
      for (int i = 0; i < n; ++i) {
         int jt = vdw::jvdw[i] - 1;
         auto iter = jmap.find(jt);
         if (iter == jmap.end()) {
            jvdwbuf[i] = jcount;
            jvec.push_back(jt);
            jmap[jt] = jcount;
            ++jcount;
         } else {
            jvdwbuf[i] = iter->second;
         }
      }

      darray::allocate(jcount * jcount, &radmin, &epsilon);

      darray::allocate(n, &vlam);

      v2scale = vdwpot::v2scale;
      v3scale = vdwpot::v3scale;
      v4scale = vdwpot::v4scale;
      v5scale = vdwpot::v5scale;

      std::vector<int> exclik;
      std::vector<real> excls;
      // see also attach.f
      const int maxn13 = 3 * sizes::maxval;
      const int maxn14 = 9 * sizes::maxval;
      const int maxn15 = 27 * sizes::maxval;
      for (int i = 0; i < n; ++i) {
         int nn;
         int bask;

         if (v2scale != 1) {
            nn = couple::n12[i];
            for (int j = 0; j < nn; ++j) {
               int k = couple::i12[i][j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(v2scale - 1);
               }
            }
         }

         if (v3scale != 1) {
            nn = couple::n13[i];
            bask = i * maxn13;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i13[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(v3scale - 1);
               }
            }
         }

         if (v4scale != 1) {
            nn = couple::n14[i];
            bask = i * maxn14;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i14[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(v4scale - 1);
               }
            }
         }

         if (v5scale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i15[bask + j];
               k -= 1;
               if (k > i) {
                  exclik.push_back(i);
                  exclik.push_back(k);
                  excls.push_back(v5scale - 1);
               }
            }
         }
      }
      nvexclude_ = excls.size();
      darray::allocate(nvexclude_, &vexclude_, &vexclude_scale_);
      darray::copyin(WAIT_NEW_Q, nvexclude_, vexclude_, exclik.data());
      darray::copyin(WAIT_NEW_Q, nvexclude_, vexclude_scale_, excls.data());

      buffer_allocate(&nev, &ev, &vir_ev);
   }

   if (op & rc_init) {
      fstr_view str = vdwpot::vdwtyp;
      if (str == "LENNARD-JONES")
         vdwtyp = evdw_t::lj;
      else if (str == "BUCKINGHAM")
         vdwtyp = evdw_t::buck;
      else if (str == "MM3-HBOND")
         vdwtyp = evdw_t::mm3hb;
      else if (str == "BUFFERED-14-7")
         vdwtyp = evdw_t::hal;
      else if (str == "GAUSSIAN")
         vdwtyp = evdw_t::gauss;
      else
         assert(false);

      ghal = vdwpot::ghal;
      dhal = vdwpot::dhal;
      scexp = mutant::scexp;
      scalpha = mutant::scalpha;
      if (static_cast<int>(evdw_t::decouple) == mutant::vcouple)
         vcouple = evdw_t::decouple;
      else if (static_cast<int>(evdw_t::annihilate) == mutant::vcouple)
         vcouple = evdw_t::annihilate;

      std::vector<int> iredbuf(n);
      std::vector<double> kredbuf(n);
      for (int i = 0; i < n; ++i) {
         int jt = vdw::ired[i] - 1;
         iredbuf[i] = jt;
         kredbuf[i] = vdw::kred[i];
      }
      darray::copyin(WAIT_NEW_Q, n, ired, iredbuf.data());
      darray::copyin(WAIT_NEW_Q, n, kred, kredbuf.data());

      darray::copyin(WAIT_NEW_Q, n, jvdw, jvdwbuf.data());
      njvdw = jcount;

      // see also kvdw.f
      std::vector<double> radvec, epsvec;
      for (int it_new = 0; it_new < jcount; ++it_new) {
         int it_old = jvec[it_new];
         int base = it_old * sizes::maxclass;
         for (int jt_new = 0; jt_new < jcount; ++jt_new) {
            int jt_old = jvec[jt_new];
            int offset = base + jt_old;
            radvec.push_back(vdw::radmin[offset]);
            epsvec.push_back(vdw::epsilon[offset]);
         }
      }
      darray::copyin(WAIT_NEW_Q, jcount * jcount, radmin, radvec.data());
      darray::copyin(WAIT_NEW_Q, jcount * jcount, epsilon, epsvec.data());

      std::vector<real> vlamvec(n);
      for (int i = 0; i < n; ++i) {
         if (mutant::mut[i]) {
            vlamvec[i] = mutant::vlambda;
         } else {
            vlamvec[i] = 1;
         }
      }
      darray::copyin(WAIT_NEW_Q, n, vlam, vlamvec.data());

      // Initialize elrc and vlrc.
      const int mode_len = 6;
      char mode[mode_len + 1] = "VDW   ";
      double elrc = 0, vlrc = 0;
      if (vdwpot::use_vcorr) {
         evcorr1(mode, &elrc, &vlrc, mode_len);
         elrc_vol = elrc * volbox();
         vlrc_vol = vlrc * volbox();
      } else {
         elrc_vol = 0;
         vlrc_vol = 0;
      }
   }
}

void evdw_reduce_xyz()
{
   evdw_reduce_xyz_acc();
}

void evdw_resolve_gradient()
{
   evdw_resolve_gradient_acc();
}

void evdw_lj(int vers)
{
   evdw_lj_acc(vers);
}

void evdw_buck(int vers)
{
   evdw_buck_acc(vers);
}

void evdw_mm3hb(int vers)
{
   evdw_mm3hb_acc(vers);
}

void evdw_hal(int vers)
{
#if TINKER_CUDART
   if (vlist_version() == NBL_SPATIAL)
      evdw_hal_cu(vers);
   else
#endif
      evdw_hal_acc(vers);
}

void evdw_gauss(int vers)
{
   evdw_gauss_acc(vers);
}

void evdw(int vers)
{
   // vdw long-range correction
   // check lrc_vol != 0 for non-PBC
   // update the global variables if !calc::analyz
   // if calc::analyz, update in energy_reduce() / virial_reduce()
   if (!(vers & calc::analyz)) {
      if ((vers & calc::energy) && elrc_vol != 0) {
         esum += elrc_vol / volbox();
      }


      if ((vers & calc::virial) && vlrc_vol != 0) {
         virial_prec term = vlrc_vol / volbox();
         vir[0] += term; // xx
         vir[4] += term; // yy
         vir[8] += term; // zz
      }
   }


   if (vdwtyp == evdw_t::lj)
      evdw_lj(vers);
   else if (vdwtyp == evdw_t::buck)
      evdw_buck(vers);
   else if (vdwtyp == evdw_t::mm3hb)
      evdw_mm3hb(vers);
   else if (vdwtyp == evdw_t::hal)
      evdw_hal(vers);
   else if (vdwtyp == evdw_t::gauss)
      evdw_gauss(vers);
   else
      assert(false);
}
TINKER_NAMESPACE_END
