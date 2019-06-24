#include "gpu/decl_dataop.h"
#include "gpu/decl_potent.h"
#include "gpu/e_strbnd.h"
#include "rc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int nstrbnd;
int (*isb)[3];
real (*sbk)[2];
real stbnunit;

real* eba;
real* vir_eba;

void estrbnd_data(int op) {
  if (!use_potent(strbnd_term))
    return;

  if (op & op_dealloc) {
    check_cudart(cudaFree(isb));
    check_cudart(cudaFree(sbk));

    check_cudart(cudaFree(eba));
    check_cudart(cudaFree(vir_eba));
  }

  if (op & op_alloc) {
    const size_t rs = sizeof(real);

    nstrbnd = count_bonded_term(strbnd_term);
    int nangle = count_bonded_term(angle_term);
    check_cudart(cudaMalloc(&isb, sizeof(int) * 3 * nangle));
    check_cudart(cudaMalloc(&sbk, rs * 2 * nangle));

    check_cudart(cudaMalloc(&eba, rs));
    check_cudart(cudaMalloc(&vir_eba, rs * 9));
  }

  if (op & op_copyin) {
    int nangle = count_bonded_term(angle_term);
    std::vector<int> ibuf(3 * nangle);
    for (int i = 0; i < 3 * nangle; ++i) {
      ibuf[i] = strbnd::isb[i] - 1;
    }
    copyin_data(&isb[0][0], ibuf.data(), 3 * nangle);
    copyin_data(&sbk[0][0], strbnd::sbk, 2 * nangle);
    stbnunit = angpot::stbnunit;
  }
}

void estrbnd(int vers) { estrbnd_acc_impl__(vers); }
}
TINKER_NAMESPACE_END
