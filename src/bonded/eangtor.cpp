#include "ff/energy.h"
#include "ff/evalence.h"
#include "ff/potent.h"
#include "math/zero.h"
#include "tool/externfunc.h"
#include <tinker/detail/angtor.hh>
#include <tinker/detail/torpot.hh>

namespace tinker {
void eangtorData(RcOp op)
{
   if (not use(Potent::ANGTOR))
      return;

   auto rc_a = rc_flag & calc::analyz;

   if (op & RcOp::DEALLOC) {
      nangtor = 0;
      darray::deallocate(iat, kant);
      if (rc_a)
         bufferDeallocate(rc_flag, eat, vir_eat, deatx, deaty, deatz);
      eat = nullptr;
      vir_eat = nullptr;
      deatx = nullptr;
      deaty = nullptr;
      deatz = nullptr;
   }

   if (op & RcOp::ALLOC) {
      nangtor = countBondedTerm(Potent::ANGTOR);
      darray::allocate(nangtor, &iat, &kant);

      eat = eng_buf;
      vir_eat = vir_buf;
      deatx = gx;
      deaty = gy;
      deatz = gz;
      if (rc_a)
         bufferAllocate(rc_flag, &eat, &vir_eat, &deatx, &deaty, &deatz);
   }

   if (op & RcOp::INIT) {
      std::vector<int> ibuf;
      ibuf.resize(4 * nangtor);
      for (int i = 0; i < 3 * nangtor; ++i)
         ibuf[i] = angtor::iat[i] - 1;
      darray::copyin(g::q0, nangtor, iat, ibuf.data());
      darray::copyin(g::q0, nangtor, kant, angtor::kant);
      waitFor(g::q0);
      atorunit = torpot::atorunit;
   }
}

TINKER_FVOID1(acc1, cu0, eangtor, int);
void eangtor(int vers)
{
   auto rc_a = rc_flag & calc::analyz;
   auto do_e = vers & calc::energy;
   auto do_v = vers & calc::virial;
   auto do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_eat, virial_eat);
      size_t bsize = bufferSize();
      if (do_e)
         darray::zero(g::q0, bsize, eat);
      if (do_v)
         darray::zero(g::q0, bsize, vir_eat);
      if (do_g)
         darray::zero(g::q0, n, deatx, deaty, deatz);
   }

   TINKER_FCALL1(acc1, cu0, eangtor, vers);

   if (rc_a) {
      if (do_e) {
         energy_eat = energyReduce(eat);
         energy_valence += energy_eat;
      }
      if (do_v) {
         virialReduce(virial_eat, vir_eat);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eat[iv];
      }
      if (do_g)
         sumGradient(gx, gy, gz, deatx, deaty, deatz);
   }
}
}
