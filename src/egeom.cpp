#include "egeom.h"
#include "mdpq.h"
#include "potent.h"
#include <tinker/detail/restrn.hh>
#include <tinker/detail/sizes.hh>


TINKER_NAMESPACE_BEGIN
void egeom_data(rc_op op)
{
   if (!use_potent(geom_term))
      return;


   if (op & rc_dealloc) {
      darray::deallocate(igfix, gfix);

      buffer_deallocate(rc_flag, eg, degx, degy, degz, vir_eg);
   }


   if (op & rc_alloc) {
      ngfix = restrn::ngfix;
      darray::allocate(ngfix, &igfix, &gfix);

      buffer_allocate(rc_flag, &eg, &degx, &degy, &degz, &vir_eg, &energy_eg);
   }


   if (op & rc_init) {
      darray::copyin(WAIT_NEW_Q, ngfix, igfix, restrn::igfix);
      darray::copyin(WAIT_NEW_Q, ngfix, gfix, restrn::gfix);
   }
}


void egeom(int vers)
{
   egeom_acc(vers);
}
TINKER_NAMESPACE_END
