#include "array.h"
#include "gpu/e_pitors.h"
#include "md.h"
#include "util_potent.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
int npitors;
int (*ipit)[6];
real* kpit;
real ptorunit;

real* ept;
real* vir_ept;

void epitors_data(rc_op op) {
  if (!use_potent(pitors_term))
    return;

  if (op & rc_dealloc) {
    dealloc_bytes(ipit);
    dealloc_bytes(kpit);

    free_ev(ept, vir_ept);
  }

  if (op & rc_alloc) {
    int ntors = count_bonded_term(torsion_term);
    alloc_bytes(&ipit, sizeof(int) * 6 * ntors);
    alloc_bytes(&kpit, sizeof(real) * ntors);

    alloc_ev(&ept, &vir_ept);
  }

  if (op & rc_init) {
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

extern void epitors_acc_impl_(int vers);
void epitors(int vers) { epitors_acc_impl_(vers); }
TINKER_NAMESPACE_END
