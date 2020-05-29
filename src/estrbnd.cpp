#include "estrbnd.h"
#include "md.h"
#include "potent.h"
#include <tinker/detail/angpot.hh>
#include <tinker/detail/strbnd.hh>

namespace tinker {
void estrbnd_data(rc_op op)
{
   if (!use_potent(strbnd_term))
      return;

   if (op & rc_dealloc) {
      darray::deallocate(isb, sbk);

      buffer_deallocate(rc_flag, eba, vir_eba);
      buffer_deallocate(rc_flag & ~calc::analyz, debax, debay, debaz);
   }

   if (op & rc_alloc) {
      int nangle = count_bonded_term(angle_term);
      darray::allocate(nangle, &isb, &sbk);

      nstrbnd = count_bonded_term(strbnd_term);
      buffer_allocate(rc_flag, &eba, &vir_eba);
      buffer_allocate(rc_flag & ~calc::analyz, &debax, &debay, &debaz);
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


   if (rc_flag & calc::analyz) {
      if (vers & calc::energy) {
         energy_eba = energy_reduce(eba);
         energy_valence += energy_eba;
      }
      if (vers & calc::virial) {
         virial_reduce(virial_eba, vir_eba);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eba[iv];
      }
   }
   if (vers & calc::analyz)
      if (vers & calc::grad)
         sum_gradient(gx, gy, gz, debax, debay, debaz);
}
}
