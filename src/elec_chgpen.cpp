#include "elec.h"
#include "empole.h"
#include "energy.h"
#include "epolar.h"
#include "glob.chglj.h"
#include "md.h"
#include "mod.vdwpot.h"
#include "nblist.h"
#include "pmestuf.h"
#include "potent.h"
#include "switch.h"
#include "tool/io_fort_str.h"
#include <tinker/detail/atoms.hh>
#include <tinker/detail/chgpot.hh>
#include <tinker/detail/couple.hh>
#include <tinker/detail/kchrge.hh>
#include <tinker/detail/limits.hh>
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/mpole.hh>


namespace tinker {
bool use_ewald()
{
   bool flag = use_energi_elec();
   return flag && limits::use_ewald;
}


//====================================================================//


void pchg_data(rc_op op)
{
   if (!use_potent(charge_term))
      return;


   if (op & rc_dealloc) {
      darray::deallocate(pchg);
   }


   if (op & rc_alloc) {
      darray::allocate(n, &pchg);
   }


   if (op & rc_init) {
      std::vector<real> pchgbuf(n);
      for (int i = 0; i < n; ++i) {
         int itype = atoms::type[i] - 1;
         pchgbuf[i] = kchrge::chg[itype];
      }
      darray::copyin(WAIT_NEW_Q, n, pchg, pchgbuf.data());
   }
}


//====================================================================//


void pole_data(rc_op op)
{
   if (!use_potent(mpole_term) && !use_potent(polar_term) &&
       !use_potent(repuls_term))
      return;


   if (op & rc_dealloc) {
      darray::deallocate(zaxis, pole, rpole, udir, uind);
      darray::deallocate(trqx, trqy, trqz, vir_trq);
   }


   if (op & rc_alloc) {
      darray::allocate(n, &zaxis, &pole, &rpole);


      if (use_potent(polar_term)) {
         darray::allocate(n, &uind, &udir);
      } else {
         uind = nullptr;
         udir = nullptr;
      }


      if (rc_flag & calc::grad) {
         darray::allocate(n, &trqx, &trqy, &trqz);
      } else {
         trqx = nullptr;
         trqy = nullptr;
         trqz = nullptr;
      }
      if (rc_flag & calc::virial) {
         darray::allocate(buffer_size(), &vir_trq);
      } else {
         vir_trq = nullptr;
      }
   }


   if (op & rc_init) {
      // Regarding chkpole routine:
      // 1. The chiralities of the atoms are unlikely to change in MD; but
      // still possible in Monte Carlo;
      // 2. yaxis values are directly copied from Tinker, and are NOT
      // subtracted by 1 becasue of the checks in chkpole.
      static_assert(sizeof(LocalFrame) == 4 * sizeof(int), "");
      std::vector<LocalFrame> zaxisbuf(n);
      for (int i = 0; i < n; ++i) {
         zaxisbuf[i].zaxis = mpole::zaxis[i] - 1;
         zaxisbuf[i].xaxis = mpole::xaxis[i] - 1;
         zaxisbuf[i].yaxis = mpole::yaxis[i];
         fstr_view str = mpole::polaxe[i];
         int val;
         if (str == "Z-Only")
            val = pole_z_only;
         else if (str == "Z-then-X")
            val = pole_z_then_x;
         else if (str == "Bisector")
            val = pole_bisector;
         else if (str == "Z-Bisect")
            val = pole_z_bisect;
         else if (str == "3-Fold")
            val = pole_3_fold;
         else
            val = pole_none;
         zaxisbuf[i].polaxe = val;
      }
      darray::copyin(WAIT_NEW_Q, n, zaxis, zaxisbuf.data());


      std::vector<double> polebuf(mpl_total * n);
      for (int i = 0; i < n; ++i) {
         int b1 = mpl_total * i;
         int b2 = mpole::maxpole * i;
         // Tinker c = 0, dx = 1, dy = 2, dz = 3
         // Tinker qxx = 4, qxy = 5, qxz = 6
         //        qyx    , qyy = 8, qyz = 9
         //        qzx    , qzy    , qzz = 12
         polebuf[b1 + mpl_pme_0] = mpole::pole[b2 + 0];
         polebuf[b1 + mpl_pme_x] = mpole::pole[b2 + 1];
         polebuf[b1 + mpl_pme_y] = mpole::pole[b2 + 2];
         polebuf[b1 + mpl_pme_z] = mpole::pole[b2 + 3];
         polebuf[b1 + mpl_pme_xx] = mpole::pole[b2 + 4];
         polebuf[b1 + mpl_pme_xy] = mpole::pole[b2 + 5];
         polebuf[b1 + mpl_pme_xz] = mpole::pole[b2 + 6];
         polebuf[b1 + mpl_pme_yy] = mpole::pole[b2 + 8];
         polebuf[b1 + mpl_pme_yz] = mpole::pole[b2 + 9];
         polebuf[b1 + mpl_pme_zz] = mpole::pole[b2 + 12];
      }
      darray::copyin(WAIT_NEW_Q, n, pole, polebuf.data());
   }
}


void mscale_data(rc_op op)
{
   if (!use_potent(mpole_term) && !use_potent(chgtrn_term))
      return;


   if (op & rc_dealloc) {
      nmexclude = 0;
      darray::deallocate(mexclude, mexclude_scale);
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
   }


   if (op & rc_init) {
   }
}

//====================================================================//


void elec_data(rc_op op)
{
   if (op & rc_init) {
      electric = chgpot::electric;
      dielec = chgpot::dielec;
   }
   rc_man pchg42{pchg_data, op};
   rc_man pole42{pole_data, op};
   rc_man mscale42{mscale_data, op};
}


//====================================================================//


void mpole_init(int vers)
{
   if (vers & calc::grad)
      darray::zero(PROCEED_NEW_Q, n, trqx, trqy, trqz);
   if (vers & calc::virial)
      darray::zero(PROCEED_NEW_Q, buffer_size(), vir_trq);


   chkpole();
   rotpole();


   if (use_ewald()) {
      rpole_to_cmp();
      if (vir_m)
         darray::zero(PROCEED_NEW_Q, buffer_size(), vir_m);
      if (pltfm_config & CU_PLTFM) {
         bool precompute_theta = (!TINKER_CU_THETA_ON_THE_FLY_GRID_MPOLE) ||
            (!TINKER_CU_THETA_ON_THE_FLY_GRID_UIND);
         if (epme_unit.valid()) {
            if (precompute_theta)
               bspline_fill(epme_unit, 3);
         }
         if (ppme_unit.valid() && (ppme_unit != epme_unit)) {
            if (precompute_theta)
               bspline_fill(ppme_unit, 2);
         }
         if (pvpme_unit.valid()) {
            if (precompute_theta)
               bspline_fill(pvpme_unit, 2);
         }
      }
   }
}


void chkpole()
{
   chkpole_acc();
}


void rotpole()
{
   rotpole_acc();
}


void torque(int vers, grad_prec* dx, grad_prec* dy, grad_prec* dz)
{
   // #if TINKER_CUDART
   //    if (pltfm_config & CU_PLTFM)
   //    else
   // #endif
   torque_acc(vers, dx, dy, dz);
}


bool amoeba_emplar(int vers)
{
   if (mplpot::use_chgpen)
      return false;
   if (rc_flag & calc::analyz)
      return false;
   if (vers & calc::analyz)
      return false;


   return use_potent(mpole_term) && use_potent(polar_term) &&
      (mlist_version() & NBL_SPATIAL);
}


bool amoeba_empole(int vers)
{
   if (mplpot::use_chgpen)
      return false;


   if (amoeba_emplar(vers))
      return false;
   return use_potent(mpole_term);
}


bool amoeba_epolar(int vers)
{
   if (mplpot::use_chgpen)
      return false;


   if (amoeba_emplar(vers))
      return false;
   return use_potent(polar_term);
}

}
