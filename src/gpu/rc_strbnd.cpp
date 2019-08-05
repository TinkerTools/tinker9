#include "array.h"
#include "gpu/e_strbnd.h"
#include "md.h"
#include "util_potent.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
int nstrbnd;
int (*isb)[3];
real (*sbk)[2];
real stbnunit;

real* eba;
real* vir_eba;

void estrbnd_data(rc_op op) {
  if (!use_potent(strbnd_term))
    return;

  if (op & rc_dealloc) {
    dealloc_bytes(isb);
    dealloc_bytes(sbk);

    free_ev(eba, vir_eba);
  }

  if (op & rc_alloc) {
    const size_t rs = sizeof(real);

    int nangle = count_bonded_term(angle_term);
    alloc_bytes(&isb, sizeof(int) * 3 * nangle);
    alloc_bytes(&sbk, rs * 2 * nangle);

    alloc_ev(&eba, &vir_eba);
  }

  if (op & rc_init) {
    nstrbnd = count_bonded_term(strbnd_term);
    int nangle = count_bonded_term(angle_term);
    std::vector<int> ibuf(3 * nangle);
    for (int i = 0; i < 3 * nangle; ++i) {
      ibuf[i] = strbnd::isb[i] - 1;
    }
    copyin_array(&isb[0][0], ibuf.data(), 3 * nangle);
    copyin_array(&sbk[0][0], strbnd::sbk, 2 * nangle);
    stbnunit = angpot::stbnunit;
  }
}

extern void estrbnd_acc_impl_(int vers);
void estrbnd(int vers) { estrbnd_acc_impl_(vers); }
TINKER_NAMESPACE_END
