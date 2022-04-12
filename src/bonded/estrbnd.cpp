#include "ff/energy.h"
#include "ff/evalence.h"
#include "ff/potent.h"
#include "math/zero.h"
#include "tool/externfunc.h"
#include <tinker/detail/angpot.hh>
#include <tinker/detail/strbnd.hh>

namespace tinker {
void estrbndData(RcOp op)
{
   if (not use(Potent::STRBND))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & RcOp::DEALLOC) {
      darray::deallocate(isb, sbk);

      if (rc_a)
         bufferDeallocate(rc_flag, eba, vir_eba, debax, debay, debaz);
      eba = nullptr;
      vir_eba = nullptr;
      debax = nullptr;
      debay = nullptr;
      debaz = nullptr;
   }

   if (op & RcOp::ALLOC) {
      int nangle = countBondedTerm(Potent::ANGLE);
      darray::allocate(nangle, &isb, &sbk);
      nstrbnd = countBondedTerm(Potent::STRBND);

      eba = eng_buf;
      vir_eba = vir_buf;
      debax = gx;
      debay = gy;
      debaz = gz;
      if (rc_a)
         bufferAllocate(rc_flag, &eba, &vir_eba, &debax, &debay, &debaz);
   }

   if (op & RcOp::INIT) {
      int nangle = countBondedTerm(Potent::ANGLE);
      std::vector<int> ibuf(3 * nangle);
      for (int i = 0; i < 3 * nangle; ++i) {
         ibuf[i] = strbnd::isb[i] - 1;
      }
      darray::copyin(g::q0, nangle, isb, ibuf.data());
      darray::copyin(g::q0, nangle, sbk, strbnd::sbk);
      waitFor(g::q0);

      stbnunit = angpot::stbnunit;
   }
}

TINKER_FVOID2(cu, 0, acc, 1, estrbnd, int);
void estrbnd(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_eba, virial_eba);
      auto bsize = bufferSize();
      if (do_e)
         darray::zero(g::q0, bsize, eba);
      if (do_v)
         darray::zero(g::q0, bsize, vir_eba);
      if (do_g)
         darray::zero(g::q0, n, debax, debay, debaz);
   }

   TINKER_FCALL2(cu, 0, acc, 1, estrbnd, vers);

   if (rc_a) {
      if (do_e) {
         energy_eba = energyReduce(eba);
         energy_valence += energy_eba;
      }
      if (do_v) {
         virialReduce(virial_eba, vir_eba);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eba[iv];
      }
      if (do_g)
         sumGradient(gx, gy, gz, debax, debay, debaz);
   }
}
}
