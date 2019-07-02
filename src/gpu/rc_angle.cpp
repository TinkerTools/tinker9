#include "gpu/decl_mdstate.h"
#include "gpu/decl_potent.h"
#include "gpu/e_angle.h"
#include "gpu/rc.h"
#include "util/fort_str.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
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

void eangle_data(rc_t rc) {
  if (!use_potent(angle_term) && !use_potent(strbnd_term) &&
      !use_potent(opbend_term))
    return;

  if (rc & rc_dealloc) {
    check_cudart(cudaFree(iang));
    check_cudart(cudaFree(ak));
    check_cudart(cudaFree(anat));

    check_cudart(cudaFree(angtyp));

    free_ev(ea, vir_ea);
  }

  if (rc & rc_alloc) {
    const size_t rs = sizeof(real);

    nangle = count_bonded_term(angle_term);
    check_cudart(cudaMalloc(&iang, sizeof(int) * nangle * 4));
    check_cudart(cudaMalloc(&ak, rs * nangle));
    check_cudart(cudaMalloc(&anat, rs * nangle));

    check_cudart(cudaMalloc(&angtyp, sizeof(int) * nangle));

    alloc_ev(&ea, &vir_ea);
  }

  if (rc & rc_copyin) {
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

extern void eangle_acc_impl__(int vers);;
void eangle(int vers) { eangle_acc_impl__(vers); }
}
TINKER_NAMESPACE_END
