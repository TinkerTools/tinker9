#include "e_strbnd.h"

#include "ext/tinker/detail/angpot.hh"
#include "ext/tinker/detail/strbnd.hh"
#include "md.h"
#include "potent.h"

TINKER_NAMESPACE_BEGIN
void estrbnd_data(rc_op op) {
  if (!use_potent(strbnd_term))
    return;

  if (op & rc_dealloc) {
    eba_handle.dealloc();
  }

  if (op & rc_alloc) {
    int nangle = count_bonded_term(angle_term);
    isb_vec.reserve(3 * nangle);
    sbk_vec.reserve(2 * nangle);

    nstrbnd = count_bonded_term(strbnd_term);
    eba_handle.alloc(nstrbnd);
  }

  if (op & rc_init) {
    int nangle = count_bonded_term(angle_term);
    std::vector<int> ibuf(3 * nangle);
    for (int i = 0; i < 3 * nangle; ++i) {
      ibuf[i] = strbnd::isb[i] - 1;
    }
    isb_vec.copyin(ibuf.data(), 3 * nangle);
    sbk_vec.copyin(strbnd::sbk, 2 * nangle);
    stbnunit = angpot::stbnunit;
  }
}

extern void estrbnd_acc_impl_(int vers);
void estrbnd(int vers) { estrbnd_acc_impl_(vers); }
TINKER_NAMESPACE_END
