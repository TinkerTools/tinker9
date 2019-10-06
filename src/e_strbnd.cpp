#include "e_strbnd.h"
#include "md.h"
#include "potent.h"
#include <ext/tinker/detail/angpot.hh>
#include <ext/tinker/detail/strbnd.hh>

TINKER_NAMESPACE_BEGIN
void estrbnd_data(rc_op op) {
  if (!use_potent(strbnd_term))
    return;

  if (op & rc_dealloc) {
    device_array::deallocate(isb, sbk);

    eba_handle.dealloc();
  }

  if (op & rc_alloc) {
    int nangle = count_bonded_term(angle_term);
    device_array::allocate(nangle, &isb, &sbk);

    nstrbnd = count_bonded_term(strbnd_term);
    eba_handle.alloc(nstrbnd);
  }

  if (op & rc_init) {
    int nangle = count_bonded_term(angle_term);
    std::vector<int> ibuf(3 * nangle);
    for (int i = 0; i < 3 * nangle; ++i) {
      ibuf[i] = strbnd::isb[i] - 1;
    }
    device_array::copyin(nangle, isb, ibuf.data());
    device_array::copyin(nangle, sbk, strbnd::sbk);

    stbnunit = angpot::stbnunit;
  }
}

extern void estrbnd_acc_impl_(int vers);
void estrbnd(int vers) { estrbnd_acc_impl_(vers); }
TINKER_NAMESPACE_END
