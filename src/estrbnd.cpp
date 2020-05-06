#include "estrbnd.h"
#include "md.h"
#include "potent.h"
#include <tinker/detail/angpot.hh>
#include <tinker/detail/strbnd.hh>

TINKER_NAMESPACE_BEGIN
void estrbnd_data(rc_op op)
{
   if (!use_potent(strbnd_term))
      return;

   if (op & rc_dealloc) {
      darray::deallocate(isb, sbk);

      buffer_deallocate(rc_flag, eba, debax, debay, debaz, vir_eba);
   }

   if (op & rc_alloc) {
      int nangle = count_bonded_term(angle_term);
      darray::allocate(nangle, &isb, &sbk);

      nstrbnd = count_bonded_term(strbnd_term);
      buffer_allocate(rc_flag, &eba, &debax, &debay, &debaz, &vir_eba,
                      &energy_eba);
   }

   if (op & rc_init) {
      int nangle = count_bonded_term(angle_term);
      std::vector<int> ibuf(3 * nangle);
      for (int i = 0; i < 3 * nangle; ++i) {
         ibuf[i] = strbnd::isb[i] - 1;
      }
      darray::copyin(WAIT_NEW_Q, nangle, isb, ibuf.data());
      darray::copyin(WAIT_NEW_Q, nangle, sbk, strbnd::sbk);

      stbnunit = angpot::stbnunit;
   }
}


void estrbnd(int vers)
{
   estrbnd_acc(vers);
}
TINKER_NAMESPACE_END
