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
    eb_handle.dealloc();
  }

  if (op & rc_alloc) {
    const size_t rs = sizeof(real);

    nbond = count_bonded_term(bond_term);
    ibnd_vec.resize(nbond * 2);
    bl_vec.resize(nbond);
    bk_vec.resize(nbond);

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
    ibnd_vec.copyin(ibndvec.data(), nbond * 2);
    bl_vec.copyin(bndstr::bl, nbond);
    bk_vec.copyin(bndstr::bk, nbond);
  }
}

extern void ebond_acc_impl_(int vers);
void ebond(int vers) { ebond_acc_impl_(vers); }
TINKER_NAMESPACE_END
