#include "gpu/e_urey.h"
#include "mod_md.h"
#include "util_array.h"
#include "util_potent.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
// module urey
int nurey;
int (*iury)[3];
real *uk, *ul;

// module urypot
real cury, qury, ureyunit;

real* eub;
real* vir_eub;

void eurey_data(rc_op op) {
  if (!use_potent(urey_term))
    return;

  if (op & rc_dealloc) {
    dealloc_array(iury);
    dealloc_array(uk);
    dealloc_array(ul);

    free_ev(eub, vir_eub);
  }

  if (op & rc_alloc) {
    const size_t rs = sizeof(real);
    int nangle = count_bonded_term(angle_term);
    alloc_array(&iury, sizeof(int) * 3 * nangle);
    alloc_array(&uk, rs * nangle);
    alloc_array(&ul, rs * nangle);

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
