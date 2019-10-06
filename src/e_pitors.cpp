#include "e_pitors.h"
#include "md.h"
#include "potent.h"
#include <ext/tinker/detail/pitors.hh>
#include <ext/tinker/detail/torpot.hh>

TINKER_NAMESPACE_BEGIN
void epitors_data (rc_op op)
{
   if (!use_potent (pitors_term))
      return;

   if (op & rc_dealloc) {
      device_array::deallocate (ipit, kpit);

      ept_handle.dealloc ();
   }

   if (op & rc_alloc) {
      int ntors = count_bonded_term (torsion_term);
      device_array::allocate (ntors, &ipit, &kpit);

      npitors = count_bonded_term (pitors_term);
      ept_handle.alloc (npitors);
   }

   if (op & rc_init) {
      int ntors = count_bonded_term (torsion_term);
      std::vector<int> ibuf (6 * ntors);
      for (int i = 0; i < 6 * ntors; ++i)
         ibuf[i] = pitors::ipit[i] - 1;
      device_array::copyin (ntors, ipit, ibuf.data ());
      device_array::copyin (ntors, kpit, pitors::kpit);
      ptorunit = torpot::ptorunit;
   }
}

extern void epitors_acc_impl_ (int vers);
void epitors (int vers)
{
   epitors_acc_impl_ (vers);
}
TINKER_NAMESPACE_END
