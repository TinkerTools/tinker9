#include "gpu/potential.h"
#include "util/format.print.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void potential_data(int op) { e_bond_data(op); }
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_gradient1() {
  m_tinker_using_namespace;
  const char* fmt = " {:20s}{:12.6f} kcal/mol\n";
  if (gpu::use_ebond()) {
    tinker_gpu_ebond_harmonic0();
    print(stdout, fmt, "Bond0", gpu::get_ebond());
    tinker_gpu_ebond_harmonic1();
    print(stdout, fmt, "Bond1", gpu::get_ebond());
    tinker_gpu_ebond_harmonic4();
    print(stdout, fmt, "Bond4", gpu::get_ebond());
    tinker_gpu_ebond_harmonic5();
    print(stdout, fmt, "Bond5", gpu::get_ebond());
    tinker_gpu_ebond_harmonic6();
    print(stdout, fmt, "Bond6", gpu::get_ebond());
  }
}
}
