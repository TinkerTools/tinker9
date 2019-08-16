#include "e_bond.h"
#include "array.h"
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
    dealloc_bytes(ibnd);
    dealloc_bytes(bl);
    dealloc_bytes(bk);

    dealloc_ev(eb, vir_eb);
  }

  if (op & rc_alloc) {
    const size_t rs = sizeof(real);

    nbond = count_bonded_term(bond_term);
    alloc_bytes(&ibnd, sizeof(int) * nbond * 2);
    alloc_bytes(&bl, rs * nbond);
    alloc_bytes(&bk, rs * nbond);

    alloc_ev(&eb, &vir_eb);
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
    copyin_array(&ibnd[0][0], ibndvec.data(), nbond * 2);
    copyin_array(bl, bndstr::bl, nbond);
    copyin_array(bk, bndstr::bk, nbond);
  }
}

extern void ebond_acc_impl_(int vers);
void ebond(int vers) { ebond_acc_impl_(vers); }
TINKER_NAMESPACE_END
