#include "etors.h"
#include "md.h"
#include "potent.h"
#include <tinker/detail/torpot.hh>
#include <tinker/detail/tors.hh>

namespace tinker {
void etors_data(rc_op op)
{
   if (!use_potent(torsion_term))
      return;

   if (op & rc_dealloc) {
      darray::deallocate(itors, tors1, tors2, tors3, tors4, tors5, tors6);

      buffer_deallocate(rc_flag, et, vir_et);
      buffer_deallocate(rc_flag & ~calc::analyz, detx, dety, detz);
   }

   if (op & rc_alloc) {
      ntors = count_bonded_term(torsion_term);
      darray::allocate(ntors, &itors, &tors1, &tors2, &tors3, &tors4, &tors5,
                       &tors6);

      buffer_allocate(rc_flag, &et, &vir_et);
      buffer_allocate(rc_flag & ~calc::analyz, &detx, &dety, &detz);
   }

   if (op & rc_init) {
      std::vector<int> ibuf(4 * ntors);
      for (int i = 0; i < 4 * ntors; ++i) {
         ibuf[i] = tors::itors[i] - 1;
      }
      darray::copyin(WAIT_NEW_Q, ntors, itors, ibuf.data());
      darray::copyin(WAIT_NEW_Q, ntors, tors1, tors::tors1);
      darray::copyin(WAIT_NEW_Q, ntors, tors2, tors::tors2);
      darray::copyin(WAIT_NEW_Q, ntors, tors3, tors::tors3);
      darray::copyin(WAIT_NEW_Q, ntors, tors4, tors::tors4);
      darray::copyin(WAIT_NEW_Q, ntors, tors5, tors::tors5);
      darray::copyin(WAIT_NEW_Q, ntors, tors6, tors::tors6);
      torsunit = torpot::torsunit;
   }
}

void etors(int vers)
{
   etors_acc(vers);


   if (rc_flag & calc::analyz) {
      if (vers & calc::energy) {
         energy_et = energy_reduce(et);
         energy_valence += energy_et;
      }
      if (vers & calc::virial) {
         virial_reduce(virial_et, vir_et);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_et[iv];
      }
   }
   if (vers & calc::analyz)
      if (vers & calc::grad)
         sum_gradient(gx, gy, gz, detx, dety, detz);
}
}
