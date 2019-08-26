#include "e_opbend.h"
#include "array.h"
#include "ext/tinker/detail/angpot.hh"
#include "ext/tinker/detail/opbend.hh"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"

TINKER_NAMESPACE_BEGIN
void eopbend_data(rc_op op) {
  if (!use_potent(opbend_term))
    return;

  if (op & rc_dealloc) {
    dealloc_bytes(iopb);
    dealloc_bytes(opbk);

    eopb_handle.dealloc();
  }

  if (op & rc_alloc) {
    int nangle = count_bonded_term(angle_term);
    alloc_bytes(&iopb, sizeof(int) * nangle);
    alloc_bytes(&opbk, sizeof(real) * nangle);

    nopbend = count_bonded_term(opbend_term);
    eopb_handle.alloc(nopbend);
  }

  if (op & rc_init) {
    fstr_view otyp = angpot::opbtyp;
    if (otyp == "W-D-C")
      opbtyp = eopbend_t::w_d_c;
    else if (otyp == "ALLINGER")
      opbtyp = eopbend_t::allinger;
    else
      assert(false);
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
