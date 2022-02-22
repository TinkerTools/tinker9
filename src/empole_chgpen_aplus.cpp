#include "empole_chgpen_aplus.h"
#include "cflux.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"
#include "tool/host_zero.h"
#include <tinker/detail/chgpen.hh>
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/potent.hh>
#include <tinker/detail/sizes.hh>


namespace tinker {
void empole_chgpen_aplus_data(rc_op op)
{
   if (not use_potent(mpole_term) and not use_potent(chgtrn_term))
      return;
   if (not mplpot::use_chgpen)
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      if (rc_a) {
         buffer_deallocate(rc_flag, nem);
         buffer_deallocate(rc_flag, em, vir_em, demx, demy, demz);
      }
      nem = nullptr;
      em = nullptr;
      vir_em = nullptr;
      demx = nullptr;
      demy = nullptr;
      demz = nullptr;
   }

   if (op & rc_alloc) {
      nem = nullptr;
      em = eng_buf_elec;
      vir_em = vir_buf_elec;
      demx = gx_elec;
      demy = gy_elec;
      demz = gz_elec;

      if (rc_a) {
         buffer_allocate(rc_flag, &nem);
         buffer_allocate(rc_flag, &em, &vir_em, &demx, &demy, &demz);
      }
   }

   if (op & rc_init) {
   }
}


void empole_chgpen_aplus(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_a = vers & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;
   int use_cf = potent::use_chgflx;
   int use_cfgrad = use_cf and do_g;


   host_zero(energy_em, virial_em);
   size_t bsize = buffer_size();
   if (rc_a) {
      if (do_a)
         darray::zero(g::q0, bsize, nem);
      if (do_e)
         darray::zero(g::q0, bsize, em);
      if (do_v)
         darray::zero(g::q0, bsize, vir_em);
      if (do_g)
         darray::zero(g::q0, n, demx, demy, demz);
   }


   if (use_cf)
      alterchg();
   mpole_init(vers);
   if (use_cfgrad) {
      zero_pot();
   }
   if (use_ewald())
      empole_chgpen_aplus_ewald(vers, use_cfgrad);
   else
      empole_chgpen_aplus_nonewald(vers, use_cfgrad);
   torque(vers, demx, demy, demz);
   if (use_cfgrad)
      dcflux(vers, demx, demy, demz, vir_em);
   if (do_v) {
      virial_buffer u2 = vir_trq;
      virial_prec v2[9];
      virial_reduce(v2, u2);
      for (int iv = 0; iv < 9; ++iv) {
         virial_em[iv] += v2[iv];
         virial_elec[iv] += v2[iv];
      }
   }


   if (rc_a) {
      if (do_e) {
         energy_buffer u = em;
         energy_prec e = energy_reduce(u);
         energy_em += e;
         energy_elec += e;
      }
      if (do_v) {
         virial_buffer u1 = vir_em;
         virial_prec v1[9];
         virial_reduce(v1, u1);
         for (int iv = 0; iv < 9; ++iv) {
            virial_em[iv] += v1[iv];
            virial_elec[iv] += v1[iv];
         }
      }
      if (do_g)
         sum_gradient(gx_elec, gy_elec, gz_elec, demx, demy, demz);
   }
}


void empole_chgpen_aplus_nonewald(int vers, int use_cf)
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      empole_chgpen_aplus_nonewald_cu(vers, use_cf);
   else
#endif
      empole_chgpen_aplus_nonewald_acc(vers, use_cf);
}


void empole_chgpen_aplus_ewald(int vers, int use_cf)
{
   empole_chgpen_aplus_ewald_real_self(vers, use_cf);
   empole_chgpen_aplus_ewald_recip(vers, use_cf);
}


void empole_chgpen_aplus_ewald_real_self(int vers, int use_cf)
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      empole_chgpen_aplus_ewald_real_self_cu(vers, use_cf);
   else
#endif
      empole_chgpen_aplus_ewald_real_self_acc(vers, use_cf);
}


void empole_chgpen_aplus_ewald_recip(int vers, int use_cf)
{
   empole_chgpen_aplus_ewald_recip_acc(vers, use_cf);
}
}
