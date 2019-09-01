#include "e_angle.h"

#include "ext/tinker/detail/angbnd.hh"
#include "ext/tinker/detail/angpot.hh"
#include "io_fort_str.h"
#include "md.h"
#include "potent.h"

TINKER_NAMESPACE_BEGIN
void eangle_data(rc_op op) {
  if (!use_potent(angle_term) && !use_potent(strbnd_term) &&
      !use_potent(opbend_term))
    return;

  if (op & rc_dealloc) {
    device_array::deallocate(iang, ak, anat, angtyp);

    ea_handle.dealloc();
  }

  if (op & rc_alloc) {
    nangle = count_bonded_term(angle_term);
    device_array::allocate(&iang, nangle);
    device_array::allocate(&ak, nangle);
    device_array::allocate(&anat, nangle);
    device_array::allocate(&angtyp, nangle);

    ea_handle.alloc(nangle);
  }

  if (op & rc_init) {
    std::vector<int> iangvec(nangle * 4);
    for (size_t i = 0; i < iangvec.size(); ++i) {
      iangvec[i] = angbnd::iang[i] - 1;
    }
    device_array::copyin(iang, iangvec.data(), nangle);
    device_array::copyin(ak, angbnd::ak, nangle);
    device_array::copyin(anat, angbnd::anat, nangle);

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
    device_array::copyin(angtyp, angtypvec.data(), nangle);
  }
}

extern void eangle_acc_impl_(int vers);
void eangle(int vers) { eangle_acc_impl_(vers); }
TINKER_NAMESPACE_END
