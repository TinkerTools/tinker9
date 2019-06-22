#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "gpu/e_angle.h"
#include "rc.h"
#include "util/fort_str.h"
#include <ext/tinker/tinker_mod.h>

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
int use_eangle() { return potent::use_angle; }

int count_eangle() {
  if (!use_eangle())
    return -1;

  return nangle;
}

void eangle_data(int op) {
  if (!use_eangle())
    return;

  if (op & op_dealloc) {
    check_cudart(cudaFree(iang));
    check_cudart(cudaFree(ak));
    check_cudart(cudaFree(anat));

    check_cudart(cudaFree(angtyp));

    check_cudart(cudaFree(ea));
    check_cudart(cudaFree(vir_ea));
  }

  if (op & op_alloc) {
    const size_t rs = sizeof(real);

    nangle = angbnd::nangle;
    check_cudart(cudaMalloc(&iang, sizeof(int) * nangle * 4));
    check_cudart(cudaMalloc(&ak, rs * nangle));
    check_cudart(cudaMalloc(&anat, rs * nangle));

    check_cudart(cudaMalloc(&angtyp, sizeof(int) * nangle));

    check_cudart(cudaMalloc(&ea, rs));
    check_cudart(cudaMalloc(&vir_ea, rs * 9));
  }

  if (op & op_copyin) {
    std::vector<int> iangvec(nangle * 4);
    for (size_t i = 0; i < iangvec.size(); ++i) {
      iangvec[i] = angbnd::iang[i] - 1;
    }
    copyin_data(&iang[0][0], iangvec.data(), nangle * 4);
    copyin_data(ak, angbnd::ak, nangle);
    copyin_data(anat, angbnd::anat, nangle);

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
    copyin_data(angtyp, iangvec.data(), nangle);
  }
}

void eangle(int vers) { eangle_acc_impl__(vers); }
}
TINKER_NAMESPACE_END
