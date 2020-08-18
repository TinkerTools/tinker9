#include "epolar.h"
#include "epolar_chgpen.h"
#include "empole_chgpen.h"
#include "md.h"
#include "nblist.h"
#include "pme.h"
#include "potent.h"
#include "tool/host_zero.h"
#include "tool/io_print.h"
#include <map>
#include <tinker/detail/chgpen.hh>
#include <tinker/detail/couple.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/polgrp.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/sizes.hh>
#include <tinker/detail/units.hh>


namespace tinker {
void epolar_chgpen_data(rc_op op)
{
   if (!use_potent(polar_term))
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
   }

   if (op & rc_init) {
      udiag = polpot::udiag;

      const double polmin = 1.0e-16;
      std::vector<double> pinvbuf(n);
      for (int i = 0; i < n; ++i) {
         pinvbuf[i] = 1.0 / std::max(polar::polarity[i], polmin);
      }
      darray::copyin(WAIT_NEW_Q, n, polarity, polar::polarity);
      darray::copyin(WAIT_NEW_Q, n, polarity_inv, pinvbuf.data());
   }
}


void induce2(real (*ud)[3])
{
   induce_mutual_pcg2(ud);

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


void epolar_chgpen(int vers)
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
      epolar_chgpen_ewald(vers);
   else
      epolar_chgpen_nonewald(vers);
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


void epolar_chgpen_nonewald(int vers)
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
         epolar_chgpen_nonewald_cu(ver2, uind);
      else
#endif
         epolar_chgpen_nonewald_acc(ver2, uind);
   }
}


void epolar_chgpen_ewald(int vers)
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
      epolar_chgpen_ewald_real(ver2);
      epolar_chgpen_ewald_recip_self(ver2);
   }
}


void epolar_chgpen_ewald_real(int vers)
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      epolar_chgpen_ewald_real_cu(vers, uind);
   else
#endif
      epolar_chgpen_ewald_real_acc(vers, uind);
}


void epolar_chgpen_ewald_recip_self(int vers)
{
   epolar_chgpen_ewald_recip_self_acc(vers, uind);
}

}
