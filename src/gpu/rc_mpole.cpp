#include "gpu/decl_dataop.h"
#include "gpu/decl_pme.h"
#include "gpu/decl_switch.h"
#include "gpu/e_mpole.h"
#include "rc_cudart.h"
#include "util/fort_str.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int empole_electyp;
std::string empole_electyp_str;

real* em;
int* nem;
real* vir_em;

int use_empole() { return potent::use_mpole; }

void get_empole_type(int& typ, std::string& typ_str) {
  if (use_ewald()) {
    typ = elec_ewald;
    typ_str = "EWALD";
  } else {
    typ = elec_coulomb;
    typ_str = "COULOMB";
  }
}

void empole_data(int op) {
  if (!use_empole())
    return;

  if (op & op_dealloc) {
    check_cudart(cudaFree(em));
    check_cudart(cudaFree(nem));
    check_cudart(cudaFree(vir_em));
  }

  if (op & op_alloc) {
    const size_t rs = sizeof(real);

    check_cudart(cudaMalloc(&em, rs));
    check_cudart(cudaMalloc(&nem, sizeof(int)));
    check_cudart(cudaMalloc(&vir_em, 9 * rs));
  }

  if (op & op_copyin) {
    get_empole_type(empole_electyp, empole_electyp_str);

    if (empole_electyp == elec_coulomb)
      switch_cut_off(switch_mpole, mpole_switch_cut, mpole_switch_off);
  }
}

void empole(int vers) {
  if (empole_electyp == elec_coulomb)
    empole_coulomb(vers);
  else if (empole_electyp == elec_ewald)
    empole_ewald(vers);
}
}
TINKER_NAMESPACE_END
