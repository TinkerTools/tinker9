#include "epolar_chgpen.h"
#include "cflux.h"
#include "induce_donly.h"
#include "md.h"
#include "mod.uprior.h"
#include "nblist.h"
#include "pme.h"
#include "potent.h"
#include "tool/host_zero.h"
#include "tool/io_fort_str.h"
#include "tool/io_print.h"
#include <map>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/potent.hh>
#include <tinker/detail/units.hh>
#include <tinker/detail/uprior.hh>

namespace tinker {
void epolar_chgpen_data(rc_op op)
{
   if (not use_potent(polar_term))
      return;
   if (not mplpot::use_chgpen)
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      darray::deallocate(polarity, polarity_inv);

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
      darray::deallocate(work01_, work02_, work03_, work04_, work05_);

      if (polpred == UPred::ASPC) {
         darray::deallocate(udalt_00, udalt_01, udalt_02, udalt_03, udalt_04,
                            udalt_05, udalt_06, udalt_07, udalt_08, udalt_09,
                            udalt_10, udalt_11, udalt_12, udalt_13, udalt_14,
                            udalt_15);
      } else if (polpred == UPred::GEAR) {
         darray::deallocate(udalt_00, udalt_01, udalt_02, udalt_03, udalt_04,
                            udalt_05);
      } else if (polpred == UPred::LSQR) {
         darray::deallocate(udalt_00, udalt_01, udalt_02, udalt_03, udalt_04,
                            udalt_05, udalt_06);
         darray::deallocate(udalt_lsqr_a, udalt_lsqr_b);
      }
      polpred = UPred::NONE;
      maxualt = 0;
      nualt = 0;
   }

   if (op & rc_alloc) {
      darray::allocate(n, &polarity, &polarity_inv);

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

      darray::allocate(n, &work01_, &work02_, &work03_, &work04_, &work05_);
      if (uprior::use_pred) {
         fstr_view predstr = uprior::polpred;
         if (predstr == "ASPC") {
            polpred = UPred::ASPC;
         } else if (predstr == "GEAR") {
            polpred = UPred::GEAR;
         } else {
            polpred = UPred::LSQR;
         }
      } else {
         polpred = UPred::NONE;
      }
      maxualt = 0;
      nualt = 0;

      if (polpred == UPred::ASPC) {
         maxualt = 16;
         darray::allocate(n, &udalt_00, &udalt_01, &udalt_02, &udalt_03,
                          &udalt_04, &udalt_05, &udalt_06, &udalt_07, &udalt_08,
                          &udalt_09, &udalt_10, &udalt_11, &udalt_12, &udalt_13,
                          &udalt_14, &udalt_15);
         darray::zero(g::q0, n, udalt_00, udalt_01, udalt_02, udalt_03,
                      udalt_04, udalt_05, udalt_06, udalt_07, udalt_08,
                      udalt_09, udalt_10, udalt_11, udalt_12, udalt_13,
                      udalt_14, udalt_15);
      } else if (polpred == UPred::GEAR) {
         maxualt = 6;
         darray::allocate(n, &udalt_00, &udalt_01, &udalt_02, &udalt_03,
                          &udalt_04, &udalt_05);
         darray::zero(g::q0, n, udalt_00, udalt_01, udalt_02, udalt_03,
                      udalt_04, udalt_05);
      } else if (polpred == UPred::LSQR) {
         maxualt = 7;
         darray::allocate(n, &udalt_00, &udalt_01, &udalt_02, &udalt_03,
                          &udalt_04, &udalt_05, &udalt_06);
         int lenb = maxualt - 1;
         int lena = lenb * lenb; // lenb*(lenb+1)/2 should be plenty.
         darray::allocate(lena, &udalt_lsqr_a);
         darray::allocate(lenb, &udalt_lsqr_b);
         darray::zero(g::q0, n, udalt_00, udalt_01, udalt_02, udalt_03,
                      udalt_04, udalt_05, udalt_06);
      }
   }

   if (op & rc_init) {
      udiag = polpot::udiag;

      const double polmin = 1.0e-16;
      std::vector<double> pinvbuf(n);
      for (int i = 0; i < n; ++i) {
         pinvbuf[i] = 1.0 / std::max(polar::polarity[i], polmin);
      }
      darray::copyin(g::q0, n, polarity, polar::polarity);
      darray::copyin(g::q0, n, polarity_inv, pinvbuf.data());
      wait_for(g::q0);
   }
}


void induce2(real (*ud)[3])
{
   induce_mutual_pcg2(ud);
   ulspred_save2(ud);

   if (inform::debug && use_potent(polar_term)) {
      std::vector<double> uindbuf;
      uindbuf.resize(3 * n);
      darray::copyout(g::q0, n, uindbuf.data(), ud);
      wait_for(g::q0);
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


void epolar_chgpen(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_a = vers & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;
   int use_cf = potent::use_chgflx;
   int use_cfgrad = use_cf and do_g;

   host_zero(energy_ep, virial_ep);
   size_t bsize = buffer_size();
   if (rc_a) {
      if (do_a)
         darray::zero(g::q0, bsize, nep);
      if (do_e)
         darray::zero(g::q0, bsize, ep);
      if (do_v)
         darray::zero(g::q0, bsize, vir_ep);
      if (do_g)
         darray::zero(g::q0, n, depx, depy, depz);
   }


   if (use_cf)
      alterchg();
   mpole_init(vers);
   if (use_cfgrad) {
      zero_pot();
   }
   if (use_ewald())
      epolar_chgpen_ewald(vers, use_cfgrad);
   else
      epolar_chgpen_nonewald(vers, use_cfgrad);
   torque(vers, depx, depy, depz);
   if (use_cfgrad)
      dcflux(vers, depx, depy, depz, vir_ep);
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


void epolar_chgpen_nonewald(int vers, int use_cf)
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

   induce2(uind);
   if (edot)
      epolar0_dotprod(uind, udir);
   if (vers != calc::v0) {
#if TINKER_CUDART
      if (mlist_version() & NBL_SPATIAL)
         epolar_chgpen_nonewald_cu(ver2, use_cf, uind);
      else
#endif
         epolar_chgpen_nonewald_acc(ver2, use_cf, uind);
   }
}


void epolar_chgpen_ewald(int vers, int use_cf)
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

   induce2(uind);
   if (edot)
      epolar0_dotprod(uind, udir);
   if (vers != calc::v0) {
      epolar_chgpen_ewald_real(ver2, use_cf);
      epolar_chgpen_ewald_recip_self(ver2, use_cf);
   }
}


void epolar_chgpen_ewald_real(int vers, int use_cf)
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      epolar_chgpen_ewald_real_cu(vers, use_cf, uind);
   else
#endif
      epolar_chgpen_ewald_real_acc(vers, use_cf, uind);
}


void epolar_chgpen_ewald_recip_self(int vers, int use_cf)
{
   epolar_chgpen_ewald_recip_self_acc(vers, use_cf, uind);
}
}
