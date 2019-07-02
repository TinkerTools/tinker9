#include "gpu/decl_mdstate.h"
#include "gpu/decl_potent.h"
#include "gpu/e_opbend.h"
#include "gpu/rc.h"
#include "util/fort_str.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
eopbend_t opbtyp;

int nopbend;
int* iopb;
real* opbk;
real opbunit;
real copb, qopb, popb, sopb;

real* eopb;
real* vir_eopb;

void eopbend_data(rc_t rc) {
  if (!use_potent(opbend_term))
    return;

  if (rc & rc_dealloc) {
    check_cudart(cudaFree(iopb));
    check_cudart(cudaFree(opbk));

    free_ev(eopb, vir_eopb);
  }

  if (rc & rc_alloc) {
    int nangle = count_bonded_term(angle_term);
    check_cudart(cudaMalloc(&iopb, sizeof(int) * nangle));
    check_cudart(cudaMalloc(&opbk, sizeof(real) * nangle));

    alloc_ev(&eopb, &vir_eopb);
  }

  if (rc & rc_copyin) {
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

void eopbend(int vers) { eopbend_acc_impl__(vers); }
}
TINKER_NAMESPACE_END
