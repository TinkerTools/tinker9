#include "ff/pchg/evalence.h"
#include "ff/potent.h"
#include "md/md.h"
#include "tool/zero.h"
#include <tinker/detail/strtor.hh>
#include <tinker/detail/torpot.hh>

namespace tinker {
void estrtorData(RcOp op)
{
   if (not use_potent(strtor_term))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      nstrtor = 0;
      darray::deallocate(ist, kst);
      if (rc_a)
         buffer_deallocate(rc_flag, ebt, vir_ebt, debtx, debty, debtz);
      ebt = nullptr;
      vir_ebt = nullptr;
      debtx = nullptr;
      debty = nullptr;
      debtz = nullptr;
   }

   if (op & rc_alloc) {
      nstrtor = count_bonded_term(strtor_term);
      darray::allocate(nstrtor, &ist, &kst);

      ebt = eng_buf;
      vir_ebt = vir_buf;
      debtx = gx;
      debty = gy;
      debtz = gz;
      if (rc_a)
         buffer_allocate(rc_flag, &ebt, &vir_ebt, &debtx, &debty, &debtz);
   }

   if (op & rc_init) {
      std::vector<int> ibuf;
      ibuf.resize(4 * nstrtor);
      for (int i = 0; i < 4 * nstrtor; ++i)
         ibuf[i] = strtor::ist[i] - 1;
      darray::copyin(g::q0, nstrtor, ist, ibuf.data());
      darray::copyin(g::q0, nstrtor, kst, strtor::kst);
      wait_for(g::q0);
      storunit = torpot::storunit;
   }
}

void estrtor(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_ebt, virial_ebt);
      size_t bsize = buffer_size();
      if (do_e)
         darray::zero(g::q0, bsize, ebt);
      if (do_v)
         darray::zero(g::q0, bsize, vir_ebt);
      if (do_g)
         darray::zero(g::q0, n, debtx, debty, debtz);
   }

   estrtor_acc(vers);

   if (rc_a) {
      if (do_e) {
         energy_ebt = energy_reduce(ebt);
         energy_valence += energy_ebt;
      }
      if (do_v) {
         virial_reduce(virial_ebt, vir_ebt);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_ebt[iv];
      }
      if (do_g) {
         sum_gradient(gx, gy, gz, debtx, debty, debtz);
      }
   }
}
}
