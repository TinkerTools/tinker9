#include "ff/energy.h"
#include "ff/pchg/evalence.h"
#include "ff/potent.h"
#include "tool/zero.h"
#include <tinker/detail/urey.hh>
#include <tinker/detail/urypot.hh>

namespace tinker {
void eureyData(RcOp op)
{
   if (!usePotent(Potent::UREY))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      darray::deallocate(iury, uk, ul);

      if (rc_a)
         bufferDeallocate(rc_flag, eub, vir_eub, deubx, deuby, deubz);
      eub = nullptr;
      vir_eub = nullptr;
      deubx = nullptr;
      deuby = nullptr;
      deubz = nullptr;
   }

   if (op & rc_alloc) {
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

   if (op & rc_init) {
      int nangle = countBondedTerm(Potent::ANGLE);
      std::vector<int> ibuf(3 * nangle);
      for (int i = 0; i < 3 * nangle; ++i)
         ibuf[i] = urey::iury[i] - 1;
      darray::copyin(g::q0, nangle, iury, ibuf.data());
      darray::copyin(g::q0, nangle, uk, urey::uk);
      darray::copyin(g::q0, nangle, ul, urey::ul);
      wait_for(g::q0);

      cury = urypot::cury;
      qury = urypot::qury;
      ureyunit = urypot::ureyunit;
   }
}

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
         sum_gradient(gx, gy, gz, deubx, deuby, deubz);
   }
}
}
