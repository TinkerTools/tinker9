#include "gpu/e_bond.h"
#include "mod_md.h"
#include "util_array.h"
#include "util_io.h"
#include "util_potent.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
int bndtyp = 0;
std::string bndtyp_str;

real cbnd, qbnd, bndunit;
int nbond = 0;
int (*ibnd)[2] = NULL;
real *bl, *bk;

real* eb;
real* vir_eb;

void ebond_data(rc_t rc) {
  if (!use_potent(bond_term) && !use_potent(strbnd_term))
    return;

  if (rc & rc_dealloc) {
    check_cudart(cudaFree(ibnd));
    check_cudart(cudaFree(bl));
    check_cudart(cudaFree(bk));

    free_ev(eb, vir_eb);
  }

  if (rc & rc_alloc) {
    const size_t rs = sizeof(real);

    nbond = count_bonded_term(bond_term);
    check_cudart(cudaMalloc(&ibnd, sizeof(int) * nbond * 2));
    check_cudart(cudaMalloc(&bl, rs * nbond));
    check_cudart(cudaMalloc(&bk, rs * nbond));

    alloc_ev(&eb, &vir_eb);
  }

  if (rc & rc_copyin) {
    fstr_view btyp = bndpot::bndtyp;
    if (btyp == "HARMONIC")
      bndtyp = bond_harmonic;
    else if (btyp == "MORSE")
      bndtyp = bond_morse;
    bndtyp_str = btyp.trim();

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
