#include "gpu/e_tors.h"
#include "mod_md.h"
#include "util_array.h"
#include "util_potent.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
int ntors;
int (*itors)[4];
real (*tors1)[4];
real (*tors2)[4];
real (*tors3)[4];
real (*tors4)[4];
real (*tors5)[4];
real (*tors6)[4];
real torsunit;

real* et;
real* vir_et;

void etors_data(rc_t rc) {
  if (!use_potent(torsion_term))
    return;

  if (rc & rc_dealloc) {
    check_cudart(cudaFree(itors));
    check_cudart(cudaFree(tors1));
    check_cudart(cudaFree(tors2));
    check_cudart(cudaFree(tors3));
    check_cudart(cudaFree(tors4));
    check_cudart(cudaFree(tors5));
    check_cudart(cudaFree(tors6));

    free_ev(et, vir_et);
  }

  if (rc & rc_alloc) {
    const size_t rs = sizeof(real);

    ntors = count_bonded_term(torsion_term);
    check_cudart(cudaMalloc(&itors, sizeof(int) * 4 * ntors));
    check_cudart(cudaMalloc(&tors1, rs * 4 * ntors));
    check_cudart(cudaMalloc(&tors2, rs * 4 * ntors));
    check_cudart(cudaMalloc(&tors3, rs * 4 * ntors));
    check_cudart(cudaMalloc(&tors4, rs * 4 * ntors));
    check_cudart(cudaMalloc(&tors5, rs * 4 * ntors));
    check_cudart(cudaMalloc(&tors6, rs * 4 * ntors));

    alloc_ev(&et, &vir_et);
  }

  if (rc & rc_copyin) {
    std::vector<int> ibuf(4 * ntors);
    for (int i = 0; i < 4 * ntors; ++i) {
      ibuf[i] = tors::itors[i] - 1;
    }
    copyin_array(&itors[0][0], ibuf.data(), 4 * ntors);
    copyin_array(&tors1[0][0], tors::tors1, 4 * ntors);
    copyin_array(&tors2[0][0], tors::tors2, 4 * ntors);
    copyin_array(&tors3[0][0], tors::tors3, 4 * ntors);
    copyin_array(&tors4[0][0], tors::tors4, 4 * ntors);
    copyin_array(&tors5[0][0], tors::tors5, 4 * ntors);
    copyin_array(&tors6[0][0], tors::tors6, 4 * ntors);
    torsunit = torpot::torsunit;
  }
}

extern void etors_acc_impl_(int vers);
void etors(int vers) { etors_acc_impl_(vers); }
TINKER_NAMESPACE_END
