#include "e_urey.h"
#include "array.h"
#include "ext/tinker/detail/urey.hh"
#include "ext/tinker/detail/urypot.hh"
#include "md.h"
#include "potent.h"

TINKER_NAMESPACE_BEGIN
void eurey_data(rc_op op) {
  if (!use_potent(urey_term))
    return;

  if (op & rc_dealloc) {
    eub_handle.dealloc();
  }

  if (op & rc_alloc) {
    int nangle = count_bonded_term(angle_term);
    iury_vec.resize(3 * nangle);
    uk_vec.resize(nangle);
    ul_vec.resize(nangle);

    nurey = count_bonded_term(urey_term);
    eub_handle.alloc(nurey);
  }

  if (op & rc_init) {
    int nangle = count_bonded_term(angle_term);
    std::vector<int> ibuf(3 * nangle);
    for (int i = 0; i < 3 * nangle; ++i)
      ibuf[i] = urey::iury[i] - 1;
    iury_vec.copyin(ibuf.data(), 3 * nangle);
    uk_vec.copyin(urey::uk, nangle);
    ul_vec.copyin(urey::ul, nangle);

    cury = urypot::cury;
    qury = urypot::qury;
    ureyunit = urypot::ureyunit;
  }
}

extern void eurey_acc_impl_(int vers);
void eurey(int vers) { eurey_acc_impl_(vers); }
TINKER_NAMESPACE_END
