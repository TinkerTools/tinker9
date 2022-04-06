#include "ff/energy.h"
#include "ff/evalence.h"
#include "ff/potent.h"
#include "math/zero.h"
#include "tool/iofortstr.h"
#include <cassert>
#include <tinker/detail/angpot.hh>
#include <tinker/detail/opbend.hh>

namespace tinker {
void eopbendData(RcOp op)
{
   if (not use(Potent::OPBEND))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & RcOp::DEALLOC) {
      darray::deallocate(iopb, opbk);

      if (rc_a)
         bufferDeallocate(rc_flag, eopb, vir_eopb, deopbx, deopby, deopbz);
      eopb = nullptr;
      vir_eopb = nullptr;
      deopbx = nullptr;
      deopby = nullptr;
      deopbz = nullptr;
   }

   if (op & RcOp::ALLOC) {
      int nangle = countBondedTerm(Potent::ANGLE);
      darray::allocate(nangle, &iopb, &opbk);
      nopbend = countBondedTerm(Potent::OPBEND);

      eopb = eng_buf;
      vir_eopb = vir_buf;
      deopbx = gx;
      deopby = gy;
      deopbz = gz;
      if (rc_a)
         bufferAllocate(rc_flag, &eopb, &vir_eopb, &deopbx, &deopby, &deopbz);
   }

   if (op & RcOp::INIT) {
      FstrView otyp = angpot::opbtyp;
      if (otyp == "W-D-C")
         opbtyp = OPBend::WDC;
      else if (otyp == "ALLINGER")
         opbtyp = OPBend::ALLINGER;
      else
         assert(false);
      int nangle = countBondedTerm(Potent::ANGLE);
      std::vector<int> ibuf(nangle);
      for (int i = 0; i < nangle; ++i)
         ibuf[i] = opbend::iopb[i] - 1;
      darray::copyin(g::q0, nangle, iopb, ibuf.data());
      darray::copyin(g::q0, nangle, opbk, opbend::opbk);
      waitFor(g::q0);
      opbunit = angpot::opbunit;
      copb = angpot::copb;
      qopb = angpot::qopb;
      popb = angpot::popb;
      sopb = angpot::sopb;
   }
}

extern void eopbend_acc(int);
void eopbend(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_eopb, virial_eopb);
      auto bsize = bufferSize();
      if (do_e)
         darray::zero(g::q0, bsize, eopb);
      if (do_v)
         darray::zero(g::q0, bsize, vir_eopb);
      if (do_g)
         darray::zero(g::q0, n, deopbx, deopby, deopbz);
   }

   eopbend_acc(vers);

   if (rc_a) {
      if (do_e) {
         energy_eopb = energyReduce(eopb);
         energy_valence += energy_eopb;
      }
      if (do_v) {
         virialReduce(virial_eopb, vir_eopb);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eopb[iv];
      }
      if (do_g)
         sumGradient(gx, gy, gz, deopbx, deopby, deopbz);
   }
}
}
