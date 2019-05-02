#include "gpu/e.bond.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
const int ebond_harmonic = 0x001;
const int ebond_morse = 0x002;

void e_bond_data(int op) {}

template <int USE, int BNDTYP>
void ebond_tmpl() {
  if_constexpr(BNDTYP & ebond_harmonic) {}
  else if_constexpr(BNDTYP & ebond_morse) {
  }
}

extern "C" {
TINKER_BONDED_GEN(tinker_gpu_ebond_harmonic, ebond_tmpl, ebond_harmonic);

TINKER_BONDED_GEN(tinker_gpu_ebond_morse, ebond_tmpl, ebond_morse);
}
}
TINKER_NAMESPACE_END
