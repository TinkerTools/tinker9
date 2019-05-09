#include "gpu/e.bond.h"
#include "tinker.mod.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int bndtyp = 0;
std::string bndtyp_str;

real cbnd, qbnd, bndunit;
int nbond = 0;
int (*ibnd)[2] = NULL;
real *bl, *bk;

real* eb;
int use_ebond() { return potent::use_bond; }

real get_ebond() {
  if (!use_ebond())
    return 0;

  real e;
  check_cudart(cudaMemcpy(&e, eb, sizeof(real), cudaMemcpyDeviceToHost));
  return e;
}

int count_ebond() {
  if (!use_ebond())
    return -1;

  return nbond;
}

void e_bond_data(int op) {
  if (!use_ebond())
    return;

  if (op == op_destroy) {
    check_cudart(cudaFree(ibnd));
    check_cudart(cudaFree(bl));
    check_cudart(cudaFree(bk));
    check_cudart(cudaFree(eb));
  }

  if (op == op_create) {
    fstr_view btyp = bndpot::bndtyp;
    if (btyp == "HARMONIC")
      bndtyp = ebond_harmonic;
    else if (btyp == "MORSE")
      bndtyp = ebond_morse;
    bndtyp_str = btyp.trim();

    cbnd = bndpot::cbnd;
    qbnd = bndpot::qbnd;
    bndunit = bndpot::bndunit;
    nbond = bndstr::nbond;

    const size_t rs = sizeof(real);
    check_cudart(cudaMalloc(&ibnd, sizeof(int) * nbond * 2));
    check_cudart(cudaMalloc(&bl, rs * nbond));
    check_cudart(cudaMalloc(&bk, rs * nbond));
    std::vector<int> ibndvec(nbond * 2);
    for (size_t i = 0; i < ibndvec.size(); ++i) {
      ibndvec[i] = bndstr::ibnd[i] - 1;
    }
    copyin_data_1(&ibnd[0][0], ibndvec.data(), nbond * 2);
    copyin_data_1(bl, bndstr::bl, nbond);
    copyin_data_1(bk, bndstr::bk, nbond);

    check_cudart(cudaMalloc(&eb, rs));
  }
}
}
TINKER_NAMESPACE_END
