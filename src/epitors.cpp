#include "epitors.h"
#include "md.h"
#include "mod.energi.h"
#include "potent.h"
#include <tinker/detail/pitors.hh>
#include <tinker/detail/torpot.hh>

namespace tinker {
void epitors_data(rc_op op)
{
   if (!use_potent(pitors_term))
      return;

   if (op & rc_dealloc) {
      darray::deallocate(ipit, kpit);

      buffer_deallocate(rc_flag, ept, vir_ept);
      buffer_deallocate(rc_flag & ~calc::analyz, deptx, depty, deptz);
   }

   if (op & rc_alloc) {
      int ntors = count_bonded_term(torsion_term);
      darray::allocate(ntors, &ipit, &kpit);

      npitors = count_bonded_term(pitors_term);
      buffer_allocate(rc_flag, &ept, &vir_ept);
      buffer_allocate(rc_flag & ~calc::analyz, &deptx, &depty, &deptz);
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
   epitors_acc(vers);


   if (rc_flag & calc::analyz) {
      if (vers & calc::energy) {
         energy_ept = energy_reduce(ept);
         energy_valence += energy_ept;
      }
      if (vers & calc::virial) {
         virial_reduce(virial_ept, vir_ept);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_ept[iv];
      }
   }
   if (vers & calc::analyz)
      if (vers & calc::grad)
         sum_gradient(gx, gy, gz, deptx, depty, deptz);
}
}
