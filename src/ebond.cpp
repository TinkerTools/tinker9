#include "ebond.h"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"
#include <tinker/detail/bndpot.hh>
#include <tinker/detail/bndstr.hh>

namespace tinker {
void ebond_data(rc_op op)
{
   if (!use_potent(bond_term) && !use_potent(strbnd_term))
      return;

   if (op & rc_dealloc) {
      darray::deallocate(ibnd, bl, bk);

      buffer_deallocate(rc_flag, eb, debx, deby, debz, vir_eb);
   }

   if (op & rc_alloc) {
      nbond = count_bonded_term(bond_term);
      darray::allocate(nbond, &ibnd, &bl, &bk);

      buffer_allocate(rc_flag, &eb, &debx, &deby, &debz, &vir_eb, &energy_eb);
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
   ebond_acc(vers);
}
}
