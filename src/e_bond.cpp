#include "e_bond.h"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"
#include <tinker/detail/bndpot.hh>
#include <tinker/detail/bndstr.hh>

TINKER_NAMESPACE_BEGIN
void ebond_data(rc_op op)
{
   if (!use_potent(bond_term) && !use_potent(strbnd_term))
      return;

   if (op & rc_dealloc) {
      device_array::deallocate(ibnd, bl, bk);

      buffer_deallocate(eb, vir_eb);
   }

   if (op & rc_alloc) {
      nbond = count_bonded_term(bond_term);
      device_array::allocate(nbond, &ibnd, &bl, &bk);

      buffer_allocate(&eb, &vir_eb);
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
      device_array::copyin(nbond, ibnd, ibndvec.data());
      device_array::copyin(nbond, bl, bndstr::bl);
      device_array::copyin(nbond, bk, bndstr::bk);
   }
}

void ebond(int vers)
{
   extern void ebond_acc(int);
   ebond_acc(vers);
}
TINKER_NAMESPACE_END
