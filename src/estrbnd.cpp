#include "evalence.h"
#include "md.h"
#include "potent.h"
#include "tool/zero.h"
#include <tinker/detail/angpot.hh>
#include <tinker/detail/strbnd.hh>

namespace tinker {
void estrbndData(RcOp op)
{
   if (!use_potent(strbnd_term))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      darray::deallocate(isb, sbk);

      if (rc_a)
         buffer_deallocate(rc_flag, eba, vir_eba, debax, debay, debaz);
      eba = nullptr;
      vir_eba = nullptr;
      debax = nullptr;
      debay = nullptr;
      debaz = nullptr;
   }

   if (op & rc_alloc) {
      int nangle = count_bonded_term(angle_term);
      darray::allocate(nangle, &isb, &sbk);
      nstrbnd = count_bonded_term(strbnd_term);

      eba = eng_buf;
      vir_eba = vir_buf;
      debax = gx;
      debay = gy;
      debaz = gz;
      if (rc_a)
         buffer_allocate(rc_flag, &eba, &vir_eba, &debax, &debay, &debaz);
   }

   if (op & rc_init) {
      int nangle = count_bonded_term(angle_term);
      std::vector<int> ibuf(3 * nangle);
      for (int i = 0; i < 3 * nangle; ++i) {
         ibuf[i] = strbnd::isb[i] - 1;
      }
      darray::copyin(g::q0, nangle, isb, ibuf.data());
      darray::copyin(g::q0, nangle, sbk, strbnd::sbk);
      wait_for(g::q0);

      stbnunit = angpot::stbnunit;
   }
}

void estrbnd(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_eba, virial_eba);
      auto bsize = buffer_size();
      if (do_e)
         darray::zero(g::q0, bsize, eba);
      if (do_v)
         darray::zero(g::q0, bsize, vir_eba);
      if (do_g)
         darray::zero(g::q0, n, debax, debay, debaz);
   }

   estrbnd_acc(vers);

   if (rc_a) {
      if (do_e) {
         energy_eba = energy_reduce(eba);
         energy_valence += energy_eba;
      }
      if (do_v) {
         virial_reduce(virial_eba, vir_eba);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eba[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, debax, debay, debaz);
   }
}
}
