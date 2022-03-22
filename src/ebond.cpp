#include "ff/pchg/evalence.h"
#include "md.h"
#include "ff/potent.h"
#include "tool/io.h"
#include "tool/zero.h"
#include <tinker/detail/bndpot.hh>
#include <tinker/detail/bndstr.hh>
#include <tinker/detail/potent.hh>

namespace tinker {
void ebondData(RcOp op)
{
   if (not use_potent(bond_term) and not use_potent(strbnd_term) and not use_potent(strtor_term) and
      not potent::use_chgflx)
      return;

   bool rc_a = rc_flag & calc::analyz;

   if (op & rc_dealloc) {
      darray::deallocate(ibnd, bl, bk);

      if (rc_a)
         buffer_deallocate(rc_flag, eb, vir_eb, debx, deby, debz);
      eb = nullptr;
      vir_eb = nullptr;
      debx = nullptr;
      deby = nullptr;
      debz = nullptr;
   }

   if (op & rc_alloc) {
      nbond = count_bonded_term(bond_term);
      darray::allocate(nbond, &ibnd, &bl, &bk);

      eb = eng_buf;
      vir_eb = vir_buf;
      debx = gx;
      deby = gy;
      debz = gz;
      if (rc_a)
         buffer_allocate(rc_flag, &eb, &vir_eb, &debx, &deby, &debz);
   }

   if (op & rc_init) {
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
      wait_for(g::q0);
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
      auto bsize = buffer_size();
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
         energy_eb = energy_reduce(eb);
         energy_valence += energy_eb;
      }
      if (do_v) {
         virial_reduce(virial_eb, vir_eb);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eb[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, debx, deby, debz);
   }
}
}
