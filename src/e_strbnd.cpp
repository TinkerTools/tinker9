#include "e_strbnd.h"
#include "array.h"
#include "ext/tinker/detail/angpot.hh"
#include "ext/tinker/detail/strbnd.hh"
#include "md.h"
#include "potent.h"

TINKER_NAMESPACE_BEGIN
void estrbnd_data(rc_op op) {
  if (!use_potent(strbnd_term))
    return;

  if (op & rc_dealloc) {
    dealloc_bytes(isb);
    dealloc_bytes(sbk);

    eba_handle.dealloc();
  }

  if (op & rc_alloc) {
    const size_t rs = sizeof(real);

    int nangle = count_bonded_term(angle_term);
    alloc_bytes(&isb, sizeof(int) * 3 * nangle);
    alloc_bytes(&sbk, rs * 2 * nangle);

    nstrbnd = count_bonded_term(strbnd_term);
    eba_handle.alloc(nstrbnd);
  }

  if (op & rc_init) {
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
