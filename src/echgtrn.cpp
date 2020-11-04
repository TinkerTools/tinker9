#include "echgtrn.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"
#include "tool/darray.h"
#include "tool/host_zero.h"
#include <tinker/detail/chgtrn.hh>
#include <tinker/detail/ctrpot.hh>


namespace tinker {
void echgtrn_data(rc_op op)
{
   if (!use_potent(chgtrn_term))
      return;


   bool rc_a = rc_flag & calc::analyz;


   if (op & rc_dealloc) {
      darray::deallocate(chgct, dmpct);

      if (rc_a) {
         buffer_deallocate(calc::analyz, nct);
         buffer_deallocate(rc_flag, ect, vir_ect, dectx, decty, dectz);
      }
      nct = nullptr;
      ect = nullptr;
      vir_ect = nullptr;
      dectx = nullptr;
      decty = nullptr;
      dectz = nullptr;
   }


   if (op & rc_alloc) {
      darray::allocate(n, &chgct, &dmpct);


      nct = nullptr;
      ect = eng_buf_elec;
      vir_ect = vir_buf_elec;
      dectx = gx_elec;
      decty = gy_elec;
      dectz = gz_elec;
      if (rc_a) {
         buffer_allocate(rc_flag, &nct);
         buffer_allocate(rc_flag, &ect, &vir_ect, &dectx, &decty, &dectz);
      }
   }


   if (op & rc_init) {
      darray::copyin(g::q0, n, chgct, chgtrn::chgct);
      std::vector<real> dmpctvec(n);
      for (int i = 0; i < n; ++i) {
         real idmp = chgtrn::dmpct[i];
         if (idmp == 0)
            idmp = 1000;
         dmpctvec[i] = idmp;
      }
      darray::copyin(g::q0, n, dmpct, dmpctvec.data());
      wait_for(g::q0);
   }
}


void echgtrn(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_a = vers & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   host_zero(energy_ect, virial_ect);
   size_t bsize = buffer_size();
   if (rc_a) {
      if (do_a)
         darray::zero(g::q0, bsize, nct);
      if (do_e)
         darray::zero(g::q0, bsize, ect);
      if (do_v)
         darray::zero(g::q0, bsize, vir_ect);
      if (do_g)
         darray::zero(g::q0, n, dectx, decty, dectz);
   }


#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      echgtrn_cu(vers);
   else
      ;
#endif


   if (rc_a) {
      if (do_e) {
         energy_buffer u = ect;
         energy_prec e = energy_reduce(u);
         energy_ect += e;
         energy_elec += e;
      }
      if (do_v) {
         virial_buffer u = vir_ect;
         virial_prec v[9];
         virial_reduce(v, u);
         for (int iv = 0; iv < 9; ++iv) {
            virial_ect[iv] += v[iv];
            virial_elec[iv] += v[iv];
         }
      }
      if (do_g)
         sum_gradient(gx_elec, gy_elec, gz_elec, dectx, decty, dectz);
   }
}
}
