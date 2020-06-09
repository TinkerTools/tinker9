#include "echgtrn.h"
#include "mdpq.h"
#include "potent.h"
#include "tool/darray.h"
#include <tinker/detail/chgtrn.hh>
#include <tinker/detail/ctrpot.hh>


namespace tinker {
void echgtrn_data(rc_op op)
{
   if (!use_potent(chgtrn_term))
      return;


   if (op & rc_dealloc) {
      darray::deallocate(chgct, dmpct);


      buffer_deallocate(rc_flag, ect, vir_ect);
      buffer_deallocate(rc_flag, dectx, decty, dectz);
   }


   if (op & rc_alloc) {
      darray::allocate(n, &chgct, &dmpct);


      buffer_allocate(rc_flag, &ect, &vir_ect);
      buffer_allocate(rc_flag, &dectx, &decty, &dectz);
   }


   if (op & rc_init) {
      darray::copyin(PROCEED_NEW_Q, n, chgct, chgtrn::chgct);
      darray::copyin(PROCEED_NEW_Q, n, dmpct, chgtrn::dmpct);
   }
}


void echgtrn(int vers) {}
}
