#include "e_tors.h"
#include "md.h"
#include "potent.h"
#include <ext/tinker/detail/torpot.hh>
#include <ext/tinker/detail/tors.hh>

TINKER_NAMESPACE_BEGIN
void etors_data (rc_op op)
{
   if (!use_potent (torsion_term))
      return;

   if (op & rc_dealloc) {
      device_array::deallocate (itors, tors1, tors2, tors3, tors4, tors5,
                                tors6);

      et_handle.dealloc ();
   }

   if (op & rc_alloc) {
      ntors = count_bonded_term (torsion_term);
      device_array::allocate (ntors, &itors, &tors1, &tors2, &tors3, &tors4,
                              &tors5, &tors6);

      et_handle.alloc (ntors);
   }

   if (op & rc_init) {
      std::vector<int> ibuf (4 * ntors);
      for (int i = 0; i < 4 * ntors; ++i) {
         ibuf[i] = tors::itors[i] - 1;
      }
      device_array::copyin (ntors, itors, ibuf.data ());
      device_array::copyin (ntors, tors1, tors::tors1);
      device_array::copyin (ntors, tors2, tors::tors2);
      device_array::copyin (ntors, tors3, tors::tors3);
      device_array::copyin (ntors, tors4, tors::tors4);
      device_array::copyin (ntors, tors5, tors::tors5);
      device_array::copyin (ntors, tors6, tors::tors6);
      torsunit = torpot::torsunit;
   }
}

extern void etors_acc_impl_ (int vers);
void etors (int vers)
{
   etors_acc_impl_ (vers);
}
TINKER_NAMESPACE_END
