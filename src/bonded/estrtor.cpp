#include "ff/energy.h"
#include "ff/evalence.h"
#include "ff/potent.h"
#include "math/zero.h"
#include <tinker/detail/strtor.hh>
#include <tinker/detail/torpot.hh>

namespace tinker {
void estrtorData(RcOp op)
{
   if (not usePotent(Potent::STRTOR))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & RcOp::DEALLOC) {
      nstrtor = 0;
      darray::deallocate(ist, kst);
      if (rc_a)
         bufferDeallocate(rc_flag, ebt, vir_ebt, debtx, debty, debtz);
      ebt = nullptr;
      vir_ebt = nullptr;
      debtx = nullptr;
      debty = nullptr;
      debtz = nullptr;
   }

   if (op & RcOp::ALLOC) {
      nstrtor = countBondedTerm(Potent::STRTOR);
      darray::allocate(nstrtor, &ist, &kst);

      ebt = eng_buf;
      vir_ebt = vir_buf;
      debtx = gx;
      debty = gy;
      debtz = gz;
      if (rc_a)
         bufferAllocate(rc_flag, &ebt, &vir_ebt, &debtx, &debty, &debtz);
   }

   if (op & RcOp::INIT) {
      std::vector<int> ibuf;
      ibuf.resize(4 * nstrtor);
      for (int i = 0; i < 4 * nstrtor; ++i)
         ibuf[i] = strtor::ist[i] - 1;
      darray::copyin(g::q0, nstrtor, ist, ibuf.data());
      darray::copyin(g::q0, nstrtor, kst, strtor::kst);
      waitFor(g::q0);
      storunit = torpot::storunit;
   }
}

extern void estrtor_acc(int);
void estrtor(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_ebt, virial_ebt);
      size_t bsize = bufferSize();
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
         energy_ebt = energyReduce(ebt);
         energy_valence += energy_ebt;
      }
      if (do_v) {
         virialReduce(virial_ebt, vir_ebt);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_ebt[iv];
      }
      if (do_g) {
         sumGradient(gx, gy, gz, debtx, debty, debtz);
      }
   }
}
}
