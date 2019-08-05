#include "array.h"
#include "gpu/e_angle.h"
#include "md.h"
#include "nblist.h"
#include "util_io.h"
#include "util_potent.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
// module angbnd
int nangle;
int (*iang)[4];
real *ak, *anat;

// module angpot
real angunit;
real cang, qang, pang, sang;
int* angtyp;

real* ea;
real* vir_ea;

void eangle_data(rc_op op) {
  if (!use_potent(angle_term) && !use_potent(strbnd_term) &&
      !use_potent(opbend_term))
    return;

  if (op & rc_dealloc) {
    dealloc_bytes(iang);
    dealloc_bytes(ak);
    dealloc_bytes(anat);

    dealloc_bytes(angtyp);

    free_ev(ea, vir_ea);
  }

  if (op & rc_alloc) {
    const size_t rs = sizeof(real);

    nangle = count_bonded_term(angle_term);
    alloc_bytes(&iang, sizeof(int) * nangle * 4);
    alloc_bytes(&ak, rs * nangle);
    alloc_bytes(&anat, rs * nangle);

    alloc_bytes(&angtyp, sizeof(int) * nangle);

    alloc_ev(&ea, &vir_ea);
  }

  if (op & rc_init) {
    std::vector<int> iangvec(nangle * 4);
    for (size_t i = 0; i < iangvec.size(); ++i) {
      iangvec[i] = angbnd::iang[i] - 1;
    }
    copyin_array(&iang[0][0], iangvec.data(), nangle * 4);
    copyin_array(ak, angbnd::ak, nangle);
    copyin_array(anat, angbnd::anat, nangle);

    angunit = angpot::angunit;
    cang = angpot::cang;
    qang = angpot::qang;
    pang = angpot::pang;
    sang = angpot::sang;
    for (int i = 0; i < nangle; ++i) {
      fstr_view atyp = angpot::angtyp[i];
      if (atyp == "IN-PLANE")
        iangvec[i] = angle_in_plane;
      else if (atyp == "HARMONIC")
        iangvec[i] = angle_harmonic;
      else if (atyp == "LINEAR")
        iangvec[i] = angle_linear;
      else if (atyp == "FOURIER")
        iangvec[i] = angle_fourier;
      else {
        assert(false);
      }
    }
    copyin_array(angtyp, iangvec.data(), nangle);
  }
}

extern void eangle_acc_impl_(int vers);
void eangle(int vers) { eangle_acc_impl_(vers); }
TINKER_NAMESPACE_END
