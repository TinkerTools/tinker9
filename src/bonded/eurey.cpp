#include "ff/energy.h"
#include "ff/evalence.h"
#include "ff/potent.h"
#include "math/zero.h"
#include <tinker/detail/urey.hh>
#include <tinker/detail/urypot.hh>

namespace tinker {
void eureyData(RcOp op)
{
   if (!use(Potent::UREY))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & RcOp::DEALLOC) {
      darray::deallocate(iury, uk, ul);

      if (rc_a)
         bufferDeallocate(rc_flag, eub, vir_eub, deubx, deuby, deubz);
      eub = nullptr;
      vir_eub = nullptr;
      deubx = nullptr;
      deuby = nullptr;
      deubz = nullptr;
   }

   if (op & RcOp::ALLOC) {
      int nangle = countBondedTerm(Potent::ANGLE);
      darray::allocate(nangle, &iury, &uk, &ul);

      nurey = countBondedTerm(Potent::UREY);

      eub = eng_buf;
      vir_eub = vir_buf;
      deubx = gx;
      deuby = gy;
      deubz = gz;
      if (rc_a)
         bufferAllocate(rc_flag, &eub, &vir_eub, &deubx, &deuby, &deubz);
   }

   if (op & RcOp::INIT) {
      int nangle = countBondedTerm(Potent::ANGLE);
      std::vector<int> ibuf(3 * nangle);
      for (int i = 0; i < 3 * nangle; ++i)
         ibuf[i] = urey::iury[i] - 1;
      darray::copyin(g::q0, nangle, iury, ibuf.data());
      darray::copyin(g::q0, nangle, uk, urey::uk);
      darray::copyin(g::q0, nangle, ul, urey::ul);
      waitFor(g::q0);

      cury = urypot::cury;
      qury = urypot::qury;
      ureyunit = urypot::ureyunit;
   }
}

extern void eurey_acc(int);
void eurey(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_eub, virial_eub);
      auto bsize = bufferSize();
      if (do_e)
         darray::zero(g::q0, bsize, eub);
      if (do_v)
         darray::zero(g::q0, bsize, vir_eub);
      if (do_g)
         darray::zero(g::q0, n, deubx, deuby, deubz);
   }

   eurey_acc(vers);

   if (rc_a) {
      if (do_e) {
         energy_eub = energyReduce(eub);
         energy_valence += energy_eub;
      }
      if (do_v) {
         virialReduce(virial_eub, vir_eub);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eub[iv];
      }
      if (do_g)
         sumGradient(gx, gy, gz, deubx, deuby, deubz);
   }
}
}
