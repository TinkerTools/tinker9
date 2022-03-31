#include "ff/energy.h"
#include "ff/pchg/evalence.h"
#include "ff/potent.h"
#include "math/zero.h"
#include <tinker/detail/imptor.hh>
#include <tinker/detail/torpot.hh>

namespace tinker {
void eimptorData(RcOp op)
{
   if (!usePotent(Potent::IMPTORS))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & RcOp::DEALLOC) {
      darray::deallocate(iitors, itors1, itors2, itors3);
      if (rc_a)
         bufferDeallocate(rc_flag, eit, vir_eit, deitx, deity, deitz);
      eit = nullptr;
      vir_eit = nullptr;
      deitx = nullptr;
      deity = nullptr;
      deitz = nullptr;
   }

   if (op & RcOp::ALLOC) {
      nitors = imptor::nitors;
      darray::allocate(nitors, &iitors, &itors1, &itors2, &itors3);
      eit = eng_buf;
      vir_eit = vir_buf;
      deitx = gx;
      deity = gy;
      deitz = gz;
      if (rc_a)
         bufferAllocate(rc_flag, &eit, &vir_eit, &deitx, &deity, &deitz);
   }

   if (op & RcOp::INIT) {
      std::vector<int> ibuf(4 * nitors);
      for (int i = 0; i < 4 * nitors; ++i) {
         ibuf[i] = imptor::iitors[i] - 1;
      }
      darray::copyin(g::q0, nitors, iitors, ibuf.data());
      darray::copyin(g::q0, nitors, itors1, imptor::itors1);
      darray::copyin(g::q0, nitors, itors2, imptor::itors2);
      darray::copyin(g::q0, nitors, itors3, imptor::itors3);
      waitFor(g::q0);
      itorunit = torpot::itorunit;
   }
}

void eimptor(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_eit, virial_eit);
      size_t bsize = bufferSize();
      if (do_e)
         darray::zero(g::q0, bsize, eit);
      if (do_v)
         darray::zero(g::q0, bsize, vir_eit);
      if (do_g)
         darray::zero(g::q0, n, deitx, deity, deitz);
   }

   eimptor_acc(vers);

   if (rc_a) {
      if (do_e) {
         energy_eit = energyReduce(eit);
         energy_valence += energy_eit;
      }
      if (do_v) {
         virialReduce(virial_eit, vir_eit);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eit[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, deitx, deity, deitz);
   }
}
}
