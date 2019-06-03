#include "gpu/decl_dataop.h"
#include "gpu/decl_switch.h"
#include "gpu/e_mpole.h"
#include "rc_cudart.h"
#include "util/fort_str.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int electyp;
std::string electyp_str;

double mpole_switch_cut, mpole_switch_off;

real* em;
int* nem;
real* vir_em;

int use_empole() { return potent::use_mpole; }

void get_empole_type(int& typ, std::string& typ_str) {
  if (limits::use_ewald) {
    typ = elec_ewald;
    typ_str = "EWALD";
  } else {
    typ = elec_coulomb;
    typ_str = "COULOMB";
  }
}
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
namespace gpu {
void e_mpole_data(int op) {
  if (!use_empole())
    return;

  if (op == op_destroy) {
    check_cudart(cudaFree(em));
    check_cudart(cudaFree(nem));
    check_cudart(cudaFree(vir_em));
  }

  if (op == op_create) {
    get_empole_type(electyp, electyp_str);

    const size_t rs = sizeof(real);

    if (electyp == elec_coulomb)
      switch_cut_off(switch_mpole, mpole_switch_cut, mpole_switch_off);
    else if (electyp == elec_ewald)
      switch_cut_off(switch_ewald, mpole_switch_cut, mpole_switch_off);

    check_cudart(cudaMalloc(&em, rs));
    check_cudart(cudaMalloc(&nem, sizeof(int)));
    check_cudart(cudaMalloc(&vir_em, 9 * rs));
  }
}
}
TINKER_NAMESPACE_END
