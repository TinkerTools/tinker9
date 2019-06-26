#include "gpu/decl_potent.h"
#include "gpu/e_pitors.h"
#include "gpu/rc.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int npitors;
int (*ipit)[6];
real* kpit;
real ptorunit;

real* ept;
real* vir_ept;

void epitors_data(rc_t rc) {
  if (!use_potent(pitors_term))
    return;

  if (rc & rc_dealloc) {
    check_cudart(cudaFree(ipit));
    check_cudart(cudaFree(kpit));

    free_ev(ept, vir_ept);
  }

  if (rc & rc_alloc) {
    int ntors = count_bonded_term(torsion_term);
    check_cudart(cudaMalloc(&ipit, sizeof(int) * 6 * ntors));
    check_cudart(cudaMalloc(&kpit, sizeof(real) * ntors));

    alloc_ev(&ept, &vir_ept);
  }

  if (rc & rc_copyin) {
    npitors = count_bonded_term(pitors_term);
    int ntors = count_bonded_term(torsion_term);
    std::vector<int> ibuf(6 * ntors);
    for (int i = 0; i < 6 * ntors; ++i)
      ibuf[i] = pitors::ipit[i] - 1;
    copyin_array(&ipit[0][0], ibuf.data(), 6 * ntors);
    copyin_array(kpit, pitors::kpit, ntors);
    ptorunit = torpot::ptorunit;
  }
}

void epitors(int vers) { epitors_acc_impl__(vers); }
}
TINKER_NAMESPACE_END
