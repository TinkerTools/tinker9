#include "ff/energy.h"
#include "ff/evalence.h"
#include "ff/potent.h"
#include "math/zero.h"
#include "tool/externfunc.h"
#include <tinker/detail/pitors.hh>
#include <tinker/detail/torpot.hh>

namespace tinker {
void epitorsData(RcOp op)
{
   if (not use(Potent::PITORS))
      return;

   auto rc_a = rc_flag & calc::analyz;

   if (op & RcOp::DEALLOC) {
      darray::deallocate(ipit, kpit);

      if (rc_a)
         bufferDeallocate(rc_flag, ept, vir_ept, deptx, depty, deptz);
      ept = nullptr;
      vir_ept = nullptr;
      deptx = nullptr;
      depty = nullptr;
      deptz = nullptr;
   }

   if (op & RcOp::ALLOC) {
      int ntors = countBondedTerm(Potent::TORSION);
      darray::allocate(ntors, &ipit, &kpit);
      npitors = countBondedTerm(Potent::PITORS);

      ept = eng_buf;
      vir_ept = vir_buf;
      deptx = gx;
      depty = gy;
      deptz = gz;
      if (rc_a)
         bufferAllocate(rc_flag, &ept, &vir_ept, &deptx, &depty, &deptz);
   }

   if (op & RcOp::INIT) {
      int ntors = countBondedTerm(Potent::TORSION);
      std::vector<int> ibuf(6 * ntors);
      for (int i = 0; i < 6 * ntors; ++i)
         ibuf[i] = pitors::ipit[i] - 1;
      darray::copyin(g::q0, ntors, ipit, ibuf.data());
      darray::copyin(g::q0, ntors, kpit, pitors::kpit);
      waitFor(g::q0);
      ptorunit = torpot::ptorunit;
   }
}

TINKER_FVOID1(acc1, cu0, epitors, int);
void epitors(int vers)
{
   auto rc_a = rc_flag & calc::analyz;
   auto do_e = vers & calc::energy;
   auto do_v = vers & calc::virial;
   auto do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_ept, virial_ept);
      auto bsize = bufferSize();
      if (do_e)
         darray::zero(g::q0, bsize, ept);
      if (do_v)
         darray::zero(g::q0, bsize, vir_ept);
      if (do_g)
         darray::zero(g::q0, n, deptx, depty, deptz);
   }

   TINKER_FCALL1(acc1, cu0, epitors, vers);

   if (rc_a) {
      if (do_e) {
         energy_ept = energyReduce(ept);
         energy_valence += energy_ept;
      }
      if (do_v) {
         virialReduce(virial_ept, vir_ept);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_ept[iv];
      }
      if (do_g)
         sumGradient(gx, gy, gz, deptx, depty, deptz);
   }
}
}
