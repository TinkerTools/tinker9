#include "e_geom.h"
#include "potent.h"
#include <tinker/detail/restrn.hh>
#include <tinker/detail/sizes.hh>


TINKER_NAMESPACE_BEGIN
void egeom_data(rc_op op)
{
   if (!use_potent(geom_term))
      return;


   if (op & rc_dealloc) {
      device_array::deallocate(igfix, gfix);

      buffer_deallocate(eg, vir_eg);
   }


   if (op & rc_alloc) {
      ngfix = restrn::ngfix;
      device_array::allocate(ngfix, &igfix, &gfix);

      buffer_allocate(&eg, &vir_eg);
   }


   if (op & rc_init) {
      device_array::copyin(WAIT_NEW_Q, ngfix, igfix, restrn::igfix);
      device_array::copyin(WAIT_NEW_Q, ngfix, gfix, restrn::gfix);
   }
}


extern void egeom_acc(int vers);
void egeom(int vers)
{
   egeom_acc(vers);
}
TINKER_NAMESPACE_END
