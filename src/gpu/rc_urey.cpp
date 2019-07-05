#include "gpu/decl_mdstate.h"
#include "gpu/decl_potent.h"
#include "gpu/e_urey.h"
#include "gpu/rc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
// module urey
int nurey;
int (*iury)[3];
real *uk, *ul;

// module urypot
real cury, qury, ureyunit;

real* eub;
real* vir_eub;

void eurey_data(rc_t rc) {
  if (!use_potent(urey_term))
    return;

  if (rc & rc_dealloc) {
    check_cudart(cudaFree(iury));
    check_cudart(cudaFree(uk));
    check_cudart(cudaFree(ul));

    free_ev(eub, vir_eub);
  }

  if (rc & rc_alloc) {
    const size_t rs = sizeof(real);
    int nangle = count_bonded_term(angle_term);
    check_cudart(cudaMalloc(&iury, sizeof(int) * 3 * nangle));
    check_cudart(cudaMalloc(&uk, rs * nangle));
    check_cudart(cudaMalloc(&ul, rs * nangle));

    alloc_ev(&eub, &vir_eub);
  }

  if (rc & rc_copyin) {
    nurey = count_bonded_term(urey_term);
    int nangle = count_bonded_term(angle_term);
    std::vector<int> ibuf(3 * nangle);
    for (int i = 0; i < 3 * nangle; ++i)
      ibuf[i] = urey::iury[i] - 1;
    copyin_array(&iury[0][0], ibuf.data(), 3 * nangle);
    copyin_array(uk, urey::uk, nangle);
    copyin_array(ul, urey::ul, nangle);

    cury = urypot::cury;
    qury = urypot::qury;
    ureyunit = urypot::ureyunit;
  }
}

extern void eurey_acc_impl_(int vers);
void eurey(int vers) { eurey_acc_impl_(vers); }
}
TINKER_NAMESPACE_END
