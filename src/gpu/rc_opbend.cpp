#include "gpu/e_opbend.h"
#include "mod_md.h"
#include "util_array.h"
#include "util_io.h"
#include "util_potent.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
eopbend_t opbtyp;

int nopbend;
int* iopb;
real* opbk;
real opbunit;
real copb, qopb, popb, sopb;

real* eopb;
real* vir_eopb;

void eopbend_data(rc_op op) {
  if (!use_potent(opbend_term))
    return;

  if (op & rc_dealloc) {
    dealloc_bytes(iopb);
    dealloc_bytes(opbk);

    free_ev(eopb, vir_eopb);
  }

  if (op & rc_alloc) {
    int nangle = count_bonded_term(angle_term);
    alloc_bytes(&iopb, sizeof(int) * nangle);
    alloc_bytes(&opbk, sizeof(real) * nangle);

    alloc_ev(&eopb, &vir_eopb);
  }

  if (op & rc_init) {
    fstr_view otyp = angpot::opbtyp;
    if (otyp == "W-D-C")
      opbtyp = opbend_w_d_c;
    else if (otyp == "ALLINGER")
      opbtyp = opbend_allinger;
    else
      assert(false);
    nopbend = count_bonded_term(opbend_term);
    int nangle = count_bonded_term(angle_term);
    std::vector<int> ibuf(nangle);
    for (int i = 0; i < nangle; ++i)
      ibuf[i] = opbend::iopb[i] - 1;
    copyin_array(iopb, ibuf.data(), nangle);
    copyin_array(opbk, opbend::opbk, nangle);
    opbunit = angpot::opbunit;
    copb = angpot::copb;
    qopb = angpot::qopb;
    popb = angpot::popb;
    sopb = angpot::sopb;
  }
}

extern void eopbend_acc_impl_(int vers);
void eopbend(int vers) { eopbend_acc_impl_(vers); }
TINKER_NAMESPACE_END
