#include "gpu/e.bond.h"
#include "tinker.mod.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
real cbnd, qbnd, bndunit;
int bndtyp = 0;
int nbond = 0;
int (*ibnd)[2] = NULL;
real *bl, *bk;
#pragma acc declare device_resident(ibnd, bl, bk)

real eb = 0;
void e_bond_data(int op) {
  if (!potent::use_bond)
    return;

  if (op == op_destroy) {
    check_cudart(cudaFree(ibnd));
    check_cudart(cudaFree(bl));
    check_cudart(cudaFree(bk));
  }

  if (op == op_create) {
    cbnd = bndpot::cbnd;
    qbnd = bndpot::qbnd;
    bndunit = bndpot::bndunit;
    fstr_view btyp = bndpot::bndtyp;
    if (btyp == "HARMONIC")
      bndtyp = ebond_harmonic;
    else if (btyp == "MORSE")
      bndtyp = ebond_morse;
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
  }
}
}
TINKER_NAMESPACE_END
