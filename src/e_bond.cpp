#include "e_bond.h"

#include "ext/tinker/detail/bndpot.hh"
#include "ext/tinker/detail/bndstr.hh"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"

TINKER_NAMESPACE_BEGIN
void ebond_data(rc_op op) {
  if (!use_potent(bond_term) && !use_potent(strbnd_term))
    return;

  if (op & rc_dealloc) {
    device_array::deallocate(ibnd, bl, bk);

    eb_handle.dealloc();
  }

  if (op & rc_alloc) {
    nbond = count_bonded_term(bond_term);
    device_array::allocate(nbond, &ibnd, &bl, &bk);

    eb_handle.alloc(nbond);
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
    device_array::copyin(ibnd, ibndvec.data(), nbond);
    device_array::copyin(bl, bndstr::bl, nbond);
    device_array::copyin(bk, bndstr::bk, nbond);
  }
}

extern void ebond_acc_impl_(int vers);
void ebond(int vers) { ebond_acc_impl_(vers); }
TINKER_NAMESPACE_END
