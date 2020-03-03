#include "e_pitors.h"
#include "md.h"
#include "potent.h"
#include <tinker/detail/pitors.hh>
#include <tinker/detail/torpot.hh>

TINKER_NAMESPACE_BEGIN
void epitors_data(rc_op op)
{
   if (!use_potent(pitors_term))
      return;

   if (op & rc_dealloc) {
      darray::deallocate(ipit, kpit);

      buffer_deallocate(ept, vir_ept);
   }

   if (op & rc_alloc) {
      int ntors = count_bonded_term(torsion_term);
      darray::allocate(ntors, &ipit, &kpit);

      npitors = count_bonded_term(pitors_term);
      buffer_allocate(&ept, &vir_ept);
   }

   if (op & rc_init) {
      int ntors = count_bonded_term(torsion_term);
      std::vector<int> ibuf(6 * ntors);
      for (int i = 0; i < 6 * ntors; ++i)
         ibuf[i] = pitors::ipit[i] - 1;
      darray::copyin(WAIT_NEW_Q, ntors, ipit, ibuf.data());
      darray::copyin(WAIT_NEW_Q, ntors, kpit, pitors::kpit);
      ptorunit = torpot::ptorunit;
   }
}

void epitors(int vers)
{
   extern void epitors_acc(int);
   epitors_acc(vers);
}
TINKER_NAMESPACE_END
