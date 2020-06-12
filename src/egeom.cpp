#include "egeom.h"
#include "mdpq.h"
#include "potent.h"
#include "tool/host_zero.h"
#include <tinker/detail/restrn.hh>
#include <tinker/detail/sizes.hh>


namespace tinker {
void egeom_data(rc_op op)
{
   if (!use_potent(geom_term))
      return;


   if (op & rc_dealloc) {
      darray::deallocate(igfix, gfix);

      buffer_deallocate(rc_flag, eg, vir_eg);
      buffer_deallocate(rc_flag, degx, degy, degz);
   }


   if (op & rc_alloc) {
      ngfix = restrn::ngfix;
      darray::allocate(ngfix, &igfix, &gfix);

      buffer_allocate(rc_flag, &eg, &vir_eg);
      buffer_allocate(rc_flag, &degx, &degy, &degz);
   }


   if (op & rc_init) {
      darray::copyin(WAIT_NEW_Q, ngfix, igfix, restrn::igfix);
      darray::copyin(WAIT_NEW_Q, ngfix, gfix, restrn::gfix);
   }
}


void egeom(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   if (rc_a) {
      host_zero(energy_eg, virial_eg);
      auto bsize = buffer_size();
      if (do_e)
         darray::zero(PROCEED_NEW_Q, bsize, eg);
      if (do_v)
         darray::zero(PROCEED_NEW_Q, bsize, vir_eg);
      if (do_g)
         darray::zero(PROCEED_NEW_Q, n, degx, degy, degz);
   }


   egeom_acc(vers);


   if (rc_a) {
      if (do_e) {
         energy_eg = energy_reduce(eg);
         energy_valence += energy_eg;
      }
      if (do_v) {
         virial_reduce(virial_eg, vir_eg);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eg[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, degx, degy, degz);
   }
}
}
