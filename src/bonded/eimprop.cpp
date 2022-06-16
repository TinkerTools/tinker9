#include "ff/energy.h"
#include "ff/evalence.h"
#include "ff/potent.h"
#include "math/zero.h"
#include "tool/externfunc.h"
#include <tinker/detail/improp.hh>
#include <tinker/detail/torpot.hh>

namespace tinker {
void eimpropData(RcOp op)
{
   if (not use(Potent::IMPROP))
      return;

   auto rc_a = rc_flag & calc::analyz;

   if (op & RcOp::DEALLOC) {}

   if (op & RcOp::ALLOC) {
      niprop = improp::niprop;
      darray::allocate(niprop, &iiprop, &kprop, &vprop);
      eid = eng_buf;
      vir_eid = vir_buf;
      deidx = gx;
      deidy = gy;
      deidz = gz;
      if (rc_a)
         bufferAllocate(rc_flag, &eid, &vir_eid, &deidx, &deidy, &deidz);
   }

   if (op & RcOp::INIT) {
      std::vector<int> ibuf(4 * niprop);
      for (int i = 0; i < 4 * niprop; ++i) {
         ibuf[i] = improp::iiprop[i] - 1;
      }
      darray::copyin(g::q0, niprop, iiprop, ibuf.data());
      darray::copyin(g::q0, niprop, kprop, improp::kprop);
      darray::copyin(g::q0, niprop, vprop, improp::vprop);
      waitFor(g::q0);
      idihunit = torpot::idihunit;
   }
}

TINKER_FVOID1(acc1, cu0, eimprop, int);
void eimprop(int vers)
{
   auto rc_a = rc_flag & calc::analyz;
   auto do_e = vers & calc::energy;
   auto do_v = vers & calc::virial;
   auto do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_eid, virial_eid);
      size_t bsize = bufferSize();
      if (do_e)
         darray::zero(g::q0, bsize, eid);
      if (do_v)
         darray::zero(g::q0, bsize, vir_eid);
      if (do_g)
         darray::zero(g::q0, n, deidx, deidy, deidz);
   }

   TINKER_FCALL1(acc1, cu0, eimprop, vers);

   if (rc_a) {
      if (do_e) {
         energy_eid = energyReduce(eid);
         energy_valence += energy_eid;
      }
      if (do_v) {
         virialReduce(virial_eid, vir_eid);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eid[iv];
      }
      if (do_g)
         sumGradient(gx, gy, gz, deidx, deidy, deidz);
   }
}
}
