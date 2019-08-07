#include "e_angle.h"
#include "array.h"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"
#include <ext/tinker/detail/angbnd.hh>
#include <ext/tinker/detail/angpot.hh>

TINKER_NAMESPACE_BEGIN
void eangle_data(rc_op op) {
  if (!use_potent(angle_term) && !use_potent(strbnd_term) &&
      !use_potent(opbend_term))
    return;

  if (op & rc_dealloc) {
    dealloc_bytes(iang);
    dealloc_bytes(ak);
    dealloc_bytes(anat);

    dealloc_bytes(angtyp);

    dealloc_ev(ea, vir_ea);
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
    std::vector<eangle_t> angtypvec(nangle);
    for (int i = 0; i < nangle; ++i) {
      fstr_view atyp = angpot::angtyp[i];
      if (atyp == "IN-PLANE")
        angtypvec[i] = eangle_t::in_plane;
      else if (atyp == "HARMONIC")
        angtypvec[i] = eangle_t::harmonic;
      else if (atyp == "LINEAR")
        angtypvec[i] = eangle_t::linear;
      else if (atyp == "FOURIER")
        angtypvec[i] = eangle_t::fourier;
      else {
        assert(false);
      }
    }
    copyin_array(reinterpret_cast<int*>(angtyp),
                 reinterpret_cast<const int*>(angtypvec.data()), nangle);
  }
}

extern void eangle_acc_impl_(int vers);
void eangle(int vers) { eangle_acc_impl_(vers); }
TINKER_NAMESPACE_END
