#include "e_urey.h"
#include "array.h"
#include "md.h"
#include "potent.h"
#include <ext/tinker/detail/urey.hh>
#include <ext/tinker/detail/urypot.hh>

TINKER_NAMESPACE_BEGIN
void eurey_data(rc_op op) {
  if (!use_potent(urey_term))
    return;

  if (op & rc_dealloc) {
    dealloc_bytes(iury);
    dealloc_bytes(uk);
    dealloc_bytes(ul);

    dealloc_ev(eub, vir_eub);
  }

  if (op & rc_alloc) {
    const size_t rs = sizeof(real);
    int nangle = count_bonded_term(angle_term);
    alloc_bytes(&iury, sizeof(int) * 3 * nangle);
    alloc_bytes(&uk, rs * nangle);
    alloc_bytes(&ul, rs * nangle);

    alloc_ev(&eub, &vir_eub);
  }

  if (op & rc_init) {
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
TINKER_NAMESPACE_END
