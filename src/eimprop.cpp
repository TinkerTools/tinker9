#include "eimprop.h"
#include "md.h"
#include "potent.h"
#include "tool/host_zero.h"
#include <tinker/detail/improp.hh>
#include <tinker/detail/torpot.hh>


namespace tinker {
void eimprop_data(rc_op op)
{
   if (not use_potent(imporp_term))
      return;


   bool rc_a = rc_flag & calc::analyz;


   if (op & rc_dealloc) {
   }


   if (op & rc_alloc) {
      niprop = improp::niprop;
      darray::allocate(niprop, &iiprop, &kprop, &vprop);
      eid = eng_buf;
      vir_eid = vir_buf;
      deidx = gx;
      deidy = gy;
      deidz = gz;
      if (rc_a)
         buffer_allocate(rc_flag, &eid, &vir_eid, &deidx, &deidy, &deidz);
   }


   if (op & rc_init) {
      std::vector<int> ibuf(4 * niprop);
      for (int i = 0; i < 4 * niprop; ++i) {
         ibuf[i] = improp::iiprop[i] - 1;
      }
      darray::copyin(g::q0, niprop, iiprop, ibuf.data());
      darray::copyin(g::q0, niprop, kprop, improp::kprop);
      darray::copyin(g::q0, niprop, vprop, improp::vprop);
      wait_for(g::q0);
      idihunit = torpot::idihunit;
   }
}


void eimprop(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   if (rc_a) {
      host_zero(energy_eid, virial_eid);
      size_t bsize = buffer_size();
      if (do_e)
         darray::zero(g::q0, bsize, eid);
      if (do_v)
         darray::zero(g::q0, bsize, vir_eid);
      if (do_g)
         darray::zero(g::q0, n, deidx, deidy, deidz);
   }


   eimprop_acc(vers);


   if (rc_a) {
      if (do_e) {
         energy_eid = energy_reduce(eid);
         energy_valence += energy_eid;
      }
      if (do_v) {
         virial_reduce(virial_eid, vir_eid);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eid[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, deidx, deidy, deidz);
   }
}
}
