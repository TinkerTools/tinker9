#include "ff/energy.h"
#include "ff/pchg/evalence.h"
#include "ff/potent.h"
#include "math/zero.h"
#include "tool/io.h"
#include <tinker/detail/bndpot.hh>
#include <tinker/detail/bndstr.hh>
#include <tinker/detail/potent.hh>

namespace tinker {
void ebondData(RcOp op)
{
   if (not usePotent(Potent::BOND) and not usePotent(Potent::STRBND) and
      not usePotent(Potent::STRTOR) and not potent::use_chgflx)
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & RcOp::DEALLOC) {
      darray::deallocate(ibnd, bl, bk);

      if (rc_a)
         bufferDeallocate(rc_flag, eb, vir_eb, debx, deby, debz);
      eb = nullptr;
      vir_eb = nullptr;
      debx = nullptr;
      deby = nullptr;
      debz = nullptr;
   }

   if (op & RcOp::ALLOC) {
      nbond = countBondedTerm(Potent::BOND);
      darray::allocate(nbond, &ibnd, &bl, &bk);

      eb = eng_buf;
      vir_eb = vir_buf;
      debx = gx;
      deby = gy;
      debz = gz;
      if (rc_a)
         bufferAllocate(rc_flag, &eb, &vir_eb, &debx, &deby, &debz);
   }

   if (op & RcOp::INIT) {
      FstrView btyp = bndpot::bndtyp;
      if (btyp == "HARMONIC")
         bndtyp = ebond_t::harmonic;
      else if (btyp == "MORSE")
         bndtyp = ebond_t::morse;

      cbnd = bndpot::cbnd;
      qbnd = bndpot::qbnd;
      bndunit = bndpot::bndunit;

      std::vector<int> ibndvec(nbond * 2);
      for (size_t i = 0; i < ibndvec.size(); ++i) {
         ibndvec[i] = bndstr::ibnd[i] - 1;
      }
      darray::copyin(g::q0, nbond, ibnd, ibndvec.data());
      darray::copyin(g::q0, nbond, bl, bndstr::bl);
      darray::copyin(g::q0, nbond, bk, bndstr::bk);
      waitFor(g::q0);
   }
}

void ebond(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;

   if (rc_a) {
      zeroOnHost(energy_eb, virial_eb);
      auto bsize = bufferSize();
      if (do_e)
         darray::zero(g::q0, bsize, eb);
      if (do_v)
         darray::zero(g::q0, bsize, vir_eb);
      if (do_g)
         darray::zero(g::q0, n, debx, deby, debz);
   }

   ebond_acc(vers);

   if (rc_a) {
      if (do_e) {
         energy_eb = energyReduce(eb);
         energy_valence += energy_eb;
      }
      if (do_v) {
         virialReduce(virial_eb, vir_eb);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eb[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, debx, deby, debz);
   }
}
}
