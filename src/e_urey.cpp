#include "e_urey.h"
#include "md.h"
#include "potent.h"
#include <ext/tinker/detail/urey.hh>
#include <ext/tinker/detail/urypot.hh>

TINKER_NAMESPACE_BEGIN
void eurey_data(rc_op op)
{
   if (!use_potent(urey_term))
      return;

   if (op & rc_dealloc) {
      device_array::deallocate(iury, uk, ul);

      buffer_deallocate(eub, vir_eub);
   }

   if (op & rc_alloc) {
      int nangle = count_bonded_term(angle_term);
      device_array::allocate(nangle, &iury, &uk, &ul);

      nurey = count_bonded_term(urey_term);
      buffer_allocate(&eub, &vir_eub);
   }

   if (op & rc_init) {
      int nangle = count_bonded_term(angle_term);
      std::vector<int> ibuf(3 * nangle);
      for (int i = 0; i < 3 * nangle; ++i)
         ibuf[i] = urey::iury[i] - 1;
      device_array::copyin(nangle, iury, ibuf.data());
      device_array::copyin(nangle, uk, urey::uk);
      device_array::copyin(nangle, ul, urey::ul);

      cury = urypot::cury;
      qury = urypot::qury;
      ureyunit = urypot::ureyunit;
   }
}

extern void eurey_acc_impl_(int vers);
void eurey(int vers)
{
   eurey_acc_impl_(vers);
}
TINKER_NAMESPACE_END
