#include "egeom.h"
#include "mdpq.h"
#include "mod.energi.h"
#include "potent.h"
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
      buffer_deallocate(rc_flag & ~calc::analyz, degx, degy, degz);
   }


   if (op & rc_alloc) {
      ngfix = restrn::ngfix;
      darray::allocate(ngfix, &igfix, &gfix);

      buffer_allocate(rc_flag, &eg, &vir_eg);
      buffer_allocate(rc_flag & ~calc::analyz, &degx, &degy, &degz);
   }


   if (op & rc_init) {
      darray::copyin(WAIT_NEW_Q, ngfix, igfix, restrn::igfix);
      darray::copyin(WAIT_NEW_Q, ngfix, gfix, restrn::gfix);
   }
}


void egeom(int vers)
{
   egeom_acc(vers);


   if (rc_flag & calc::analyz) {
      if (vers & calc::energy) {
         energy_eg = energy_reduce(eg);
         energy_valence += energy_eg;
      }
      if (vers & calc::virial) {
         virial_reduce(virial_eg, vir_eg);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eg[iv];
      }
   }
   if (vers & calc::analyz)
      if (vers & calc::grad)
         sum_gradient(gx, gy, gz, degx, degy, degz);
}
}
