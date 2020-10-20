#include "ebond.h"
#include "md.h"
#include "potent.h"
#include "tool/host_zero.h"
#include "tool/io_fort_str.h"
#include <tinker/detail/bndpot.hh>
#include <tinker/detail/bndstr.hh>
#include <tinker/detail/potent.hh>

namespace tinker {
void ebond_data(rc_op op)
{
   if (!use_potent(bond_term) && !use_potent(strbnd_term) &&
       !potent::use_chgflx)
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
      fstr_view btyp = bndpot::bndtyp;
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
      darray::copyin(WAIT_NEW_Q, nbond, ibnd, ibndvec.data());
      darray::copyin(WAIT_NEW_Q, nbond, bl, bndstr::bl);
      darray::copyin(WAIT_NEW_Q, nbond, bk, bndstr::bk);
   }
}

void ebond(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   if (rc_a) {
      host_zero(energy_eb, virial_eb);
      auto bsize = buffer_size();
      if (do_e)
         darray::zero(PROCEED_NEW_Q, bsize, eb);
      if (do_v)
         darray::zero(PROCEED_NEW_Q, bsize, vir_eb);
      if (do_g)
         darray::zero(PROCEED_NEW_Q, n, debx, deby, debz);
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
