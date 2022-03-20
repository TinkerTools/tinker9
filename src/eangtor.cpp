#include "evalence.h"
#include "md.h"
#include "potent.h"
#include "tool/zero.h"
#include <tinker/detail/angtor.hh>
#include <tinker/detail/torpot.hh>

namespace tinker {
void eangtorData(RcOp op)
{
   if (not use_potent(angtor_term))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      nangtor = 0;
      darray::deallocate(iat, kant);
      if (rc_a)
         buffer_deallocate(rc_flag, eat, vir_eat, deatx, deaty, deatz);
      eat = nullptr;
      vir_eat = nullptr;
      deatx = nullptr;
      deaty = nullptr;
      deatz = nullptr;
   }

   if (op & rc_alloc) {
      nangtor = count_bonded_term(angtor_term);
      darray::allocate(nangtor, &iat, &kant);

      eat = eng_buf;
      vir_eat = vir_buf;
      deatx = gx;
      deaty = gy;
      deatz = gz;
      if (rc_a)
         buffer_allocate(rc_flag, &eat, &vir_eat, &deatx, &deaty, &deatz);
   }

   if (op & rc_init) {
      std::vector<int> ibuf;
      ibuf.resize(4 * nangtor);
      for (int i = 0; i < 3 * nangtor; ++i)
         ibuf[i] = angtor::iat[i] - 1;
      darray::copyin(g::q0, nangtor, iat, ibuf.data());
      darray::copyin(g::q0, nangtor, kant, angtor::kant);
      wait_for(g::q0);
      atorunit = torpot::atorunit;
   }
}

void eangtor(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_eat, virial_eat);
      size_t bsize = buffer_size();
      if (do_e)
         darray::zero(g::q0, bsize, eat);
      if (do_v)
         darray::zero(g::q0, bsize, vir_eat);
      if (do_g)
         darray::zero(g::q0, n, deatx, deaty, deatz);
   }

   eangtor_acc(vers);

   if (rc_a) {
      if (do_e) {
         energy_eat = energy_reduce(eat);
         energy_valence += energy_eat;
      }
      if (do_v) {
         virial_reduce(virial_eat, vir_eat);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eat[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, deatx, deaty, deatz);
   }
}
}
