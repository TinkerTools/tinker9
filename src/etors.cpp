#include "etors.h"
#include "glob.accasync.h"
#include "md.h"
#include "potent.h"
#include "tool/host_zero.h"
#include <tinker/detail/torpot.hh>
#include <tinker/detail/tors.hh>

namespace tinker {
void etors_data(rc_op op)
{
   if (not use_potent(torsion_term) and not use_potent(strtor_term) and
       not use_potent(angtor_term))
      return;


   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      darray::deallocate(itors, tors1, tors2, tors3, tors4, tors5, tors6);

      if (rc_a)
         buffer_deallocate(rc_flag, et, vir_et, detx, dety, detz);
      et = nullptr;
      vir_et = nullptr;
      detx = nullptr;
      dety = nullptr;
      detz = nullptr;
   }

   if (op & rc_alloc) {
      ntors = count_bonded_term(torsion_term);
      darray::allocate(ntors, &itors, &tors1, &tors2, &tors3, &tors4, &tors5,
                       &tors6);

      et = eng_buf;
      vir_et = vir_buf;
      detx = gx;
      dety = gy;
      detz = gz;
      if (rc_a)
         buffer_allocate(rc_flag, &et, &vir_et, &detx, &dety, &detz);
   }

   if (op & rc_init) {
      std::vector<int> ibuf(4 * ntors);
      for (int i = 0; i < 4 * ntors; ++i) {
         ibuf[i] = tors::itors[i] - 1;
      }
      darray::copyin(g::q0, ntors, itors, ibuf.data());
      darray::copyin(g::q0, ntors, tors1, tors::tors1);
      darray::copyin(g::q0, ntors, tors2, tors::tors2);
      darray::copyin(g::q0, ntors, tors3, tors::tors3);
      darray::copyin(g::q0, ntors, tors4, tors::tors4);
      darray::copyin(g::q0, ntors, tors5, tors::tors5);
      darray::copyin(g::q0, ntors, tors6, tors::tors6);
      wait_for(g::q0);
      torsunit = torpot::torsunit;
   }
}

void etors(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   if (rc_a) {
      host_zero(energy_et, virial_et);
      auto bsize = buffer_size();
      if (do_e)
         darray::zero(g::q0, bsize, et);
      if (do_v)
         darray::zero(g::q0, bsize, vir_et);
      if (do_g)
         darray::zero(g::q0, n, detx, dety, detz);
   }


   etors_acc(vers);


   if (rc_a) {
      if (do_e) {
         energy_et = energy_reduce(et);
         energy_valence += energy_et;
      }
      if (do_v) {
         virial_reduce(virial_et, vir_et);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_et[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, detx, dety, detz);
   }
}
}
