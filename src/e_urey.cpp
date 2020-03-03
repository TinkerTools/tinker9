#include "e_urey.h"
#include "md.h"
#include "potent.h"
#include <tinker/detail/urey.hh>
#include <tinker/detail/urypot.hh>

TINKER_NAMESPACE_BEGIN
void eurey_data(rc_op op)
{
   if (!use_potent(urey_term))
      return;

   if (op & rc_dealloc) {
      darray::deallocate(iury, uk, ul);

      buffer_deallocate(eub, vir_eub);
   }

   if (op & rc_alloc) {
      int nangle = count_bonded_term(angle_term);
      darray::allocate(nangle, &iury, &uk, &ul);

      nurey = count_bonded_term(urey_term);
      buffer_allocate(&eub, &vir_eub);
   }

   if (op & rc_init) {
      int nangle = count_bonded_term(angle_term);
      std::vector<int> ibuf(3 * nangle);
      for (int i = 0; i < 3 * nangle; ++i)
         ibuf[i] = urey::iury[i] - 1;
      darray::copyin(WAIT_NEW_Q, nangle, iury, ibuf.data());
      darray::copyin(WAIT_NEW_Q, nangle, uk, urey::uk);
      darray::copyin(WAIT_NEW_Q, nangle, ul, urey::ul);

      cury = urypot::cury;
      qury = urypot::qury;
      ureyunit = urypot::ureyunit;
   }
}

void eurey(int vers)
{
   extern void eurey_acc(int);
   eurey_acc(vers);
}
TINKER_NAMESPACE_END
