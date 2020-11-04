#include "epitors.h"
#include "md.h"
#include "potent.h"
#include "tool/host_zero.h"
#include <tinker/detail/pitors.hh>
#include <tinker/detail/torpot.hh>

namespace tinker {
void epitors_data(rc_op op)
{
   if (!use_potent(pitors_term))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      darray::deallocate(ipit, kpit);

      if (rc_a)
         buffer_deallocate(rc_flag, ept, vir_ept, deptx, depty, deptz);
      ept = nullptr;
      vir_ept = nullptr;
      deptx = nullptr;
      depty = nullptr;
      deptz = nullptr;
   }

   if (op & rc_alloc) {
      int ntors = count_bonded_term(torsion_term);
      darray::allocate(ntors, &ipit, &kpit);
      npitors = count_bonded_term(pitors_term);

      ept = eng_buf;
      vir_ept = vir_buf;
      deptx = gx;
      depty = gy;
      deptz = gz;
      if (rc_a)
         buffer_allocate(rc_flag, &ept, &vir_ept, &deptx, &depty, &deptz);
   }

   if (op & rc_init) {
      int ntors = count_bonded_term(torsion_term);
      std::vector<int> ibuf(6 * ntors);
      for (int i = 0; i < 6 * ntors; ++i)
         ibuf[i] = pitors::ipit[i] - 1;
      darray::copyin(g::q0, ntors, ipit, ibuf.data());
      darray::copyin(g::q0, ntors, kpit, pitors::kpit);
      wait_for(g::q0);
      ptorunit = torpot::ptorunit;
   }
}

void epitors(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   if (rc_a) {
      host_zero(energy_ept, virial_ept);
      auto bsize = buffer_size();
      if (do_e)
         darray::zero(g::q0, bsize, ept);
      if (do_v)
         darray::zero(g::q0, bsize, vir_ept);
      if (do_g)
         darray::zero(g::q0, n, deptx, depty, deptz);
   }


   epitors_acc(vers);


   if (rc_a) {
      if (do_e) {
         energy_ept = energy_reduce(ept);
         energy_valence += energy_ept;
      }
      if (do_v) {
         virial_reduce(virial_ept, vir_ept);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_ept[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, deptx, depty, deptz);
   }
}
}
