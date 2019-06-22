#include "gpu/decl_dataop.h"
#include "gpu/e_bond.h"
#include "rc.h"
#include "util/fort_str.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int bndtyp = 0;
std::string bndtyp_str;

real cbnd, qbnd, bndunit;
int nbond = 0;
int (*ibnd)[2] = NULL;
real *bl, *bk;

real* eb;
real* vir_eb;
int use_ebond() { return potent::use_bond; }

int count_ebond() {
  if (!use_ebond())
    return -1;

  return nbond;
}

void ebond_data(int op) {
  if (!use_ebond())
    return;

  if (op & op_dealloc) {
    check_cudart(cudaFree(ibnd));
    check_cudart(cudaFree(bl));
    check_cudart(cudaFree(bk));
    check_cudart(cudaFree(eb));
    check_cudart(cudaFree(vir_eb));
  }

  if (op & op_alloc) {
    const size_t rs = sizeof(real);
    check_cudart(cudaMalloc(&ibnd, sizeof(int) * nbond * 2));
    check_cudart(cudaMalloc(&bl, rs * nbond));
    check_cudart(cudaMalloc(&bk, rs * nbond));
    check_cudart(cudaMalloc(&eb, rs));
    check_cudart(cudaMalloc(&vir_eb, rs * 9));
  }

  if (op & op_copyin) {
    fstr_view btyp = bndpot::bndtyp;
    if (btyp == "HARMONIC")
      bndtyp = bond_harmonic;
    else if (btyp == "MORSE")
      bndtyp = bond_morse;
    bndtyp_str = btyp.trim();

    cbnd = bndpot::cbnd;
    qbnd = bndpot::qbnd;
    bndunit = bndpot::bndunit;
    nbond = bndstr::nbond;

    std::vector<int> ibndvec(nbond * 2);
    for (size_t i = 0; i < ibndvec.size(); ++i) {
      ibndvec[i] = bndstr::ibnd[i] - 1;
    }
    copyin_data(&ibnd[0][0], ibndvec.data(), nbond * 2);
    copyin_data(bl, bndstr::bl, nbond);
    copyin_data(bk, bndstr::bk, nbond);
  }
}

void ebond(int vers) {
  if (bndtyp == bond_harmonic)
    ebond_harmonic(vers);
  else if (bndtyp == bond_morse)
    ebond_morse(vers);
}
}
TINKER_NAMESPACE_END
