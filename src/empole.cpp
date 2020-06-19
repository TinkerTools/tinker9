#include "empole.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"
#include "tool/host_zero.h"
#include <tinker/detail/sizes.hh>


namespace tinker {
void empole_data(rc_op op)
{
   if (!use_potent(mpole_term))
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


void empole(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_a = vers & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   host_zero(energy_em, virial_em);
   size_t bsize = buffer_size();
   if (rc_a) {
      if (do_a)
         darray::zero(PROCEED_NEW_Q, bsize, nem);
      if (do_e)
         darray::zero(PROCEED_NEW_Q, bsize, em);
      if (do_v)
         darray::zero(PROCEED_NEW_Q, bsize, vir_em);
      if (do_g)
         darray::zero(PROCEED_NEW_Q, n, demx, demy, demz);
   }


   mpole_init(vers);
   if (use_ewald())
      empole_ewald(vers);
   else
      empole_nonewald(vers);
   torque(vers, demx, demy, demz);
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
}
