#include "ff/energy.h"
#include "ff/pchg/evalence.h"
#include "ff/potent.h"
#include "math/zero.h"
#include <tinker/detail/restrn.hh>
#include <tinker/detail/sizes.hh>

namespace tinker {
void egeomData(RcOp op)
{
   if (!usePotent(Potent::GEOM))
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & RcOp::DEALLOC) {
      npfix = 0;
      darray::deallocate(ipfix, kpfix, xpfix, ypfix, zpfix, pfix);
      ngfix = 0;
      darray::deallocate(igfix, gfix);
      ndfix = 0;
      darray::deallocate(idfix, dfix);
      nafix = 0;
      darray::deallocate(iafix, afix);
      ntfix = 0;
      darray::deallocate(itfix, tfix);

      if (rc_a)
         bufferDeallocate(rc_flag, eg, vir_eg, degx, degy, degz);
      eg = nullptr;
      vir_eg = nullptr;
      degx = nullptr;
      degy = nullptr;
      degz = nullptr;
   }

   if (op & RcOp::ALLOC) {
      npfix = restrn::npfix;
      darray::allocate(npfix, &ipfix, &kpfix, &xpfix, &ypfix, &zpfix, &pfix);
      ngfix = restrn::ngfix;
      darray::allocate(ngfix, &igfix, &gfix);
      ndfix = restrn::ndfix;
      darray::allocate(ndfix, &idfix, &dfix);
      nafix = restrn::nafix;
      darray::allocate(nafix, &iafix, &afix);
      ntfix = restrn::ntfix;
      darray::allocate(ntfix, &itfix, &tfix);

      eg = eng_buf;
      vir_eg = vir_buf;
      degx = gx;
      degy = gy;
      degz = gz;
      if (rc_a)
         bufferAllocate(rc_flag, &eg, &vir_eg, &degx, &degy, &degz);
   }

   if (op & RcOp::INIT) {
      std::vector<int> ipfixbuf(npfix);
      for (int i = 0; i < npfix; ++i) {
         ipfixbuf[i] = restrn::ipfix[i] - 1;
      }
      darray::copyin(g::q0, npfix, ipfix, ipfixbuf.data());
      darray::copyin(g::q0, npfix, kpfix, restrn::kpfix);
      darray::copyin(g::q0, npfix, xpfix, restrn::xpfix);
      darray::copyin(g::q0, npfix, ypfix, restrn::ypfix);
      darray::copyin(g::q0, npfix, zpfix, restrn::zpfix);
      darray::copyin(g::q0, npfix, pfix, restrn::pfix);
      darray::copyin(g::q0, ngfix, igfix, restrn::igfix);
      darray::copyin(g::q0, ngfix, gfix, restrn::gfix);
      std::vector<int> idfixbuf(2 * ndfix);
      for (int i = 0; i < 2 * ndfix; ++i) {
         idfixbuf[i] = restrn::idfix[i] - 1;
      }
      darray::copyin(g::q0, ndfix, idfix, idfixbuf.data());
      darray::copyin(g::q0, ndfix, dfix, restrn::dfix);
      std::vector<int> iafixbuf(3 * nafix);
      for (int i = 0; i < 3 * nafix; ++i) {
         iafixbuf[i] = restrn::iafix[i] - 1;
      }
      darray::copyin(g::q0, nafix, iafix, iafixbuf.data());
      darray::copyin(g::q0, nafix, afix, restrn::afix);
      std::vector<int> itfixbuf(4 * ntfix);
      for (int i = 0; i < 4 * ntfix; ++i) {
         itfixbuf[i] = restrn::itfix[i] - 1;
      }
      darray::copyin(g::q0, ntfix, itfix, itfixbuf.data());
      darray::copyin(g::q0, ntfix, tfix, restrn::tfix);
      waitFor(g::q0);
   }
}

void egeom(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_eg, virial_eg);
      auto bsize = bufferSize();
      if (do_e)
         darray::zero(g::q0, bsize, eg);
      if (do_v)
         darray::zero(g::q0, bsize, vir_eg);
      if (do_g)
         darray::zero(g::q0, n, degx, degy, degz);
   }

   egeom_acc(vers);

   if (rc_a) {
      if (do_e) {
         energy_eg = energyReduce(eg);
         energy_valence += energy_eg;
      }
      if (do_v) {
         virialReduce(virial_eg, vir_eg);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eg[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, degx, degy, degz);
   }
}
}
