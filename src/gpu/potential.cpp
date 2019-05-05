#include "gpu/potential.h"
#include "gpu/acc.h"
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

  gpu::async_launches_begin(&gpu::queue_b);

  if (gpu::use_ebond()) {
    tinker_gpu_ebond_harmonic0();
    tinker_gpu_ebond_harmonic4();
    tinker_gpu_ebond_harmonic5();
    tinker_gpu_ebond_harmonic6();
    tinker_gpu_ebond_harmonic1();
  }

  gpu::async_launches_end(gpu::queue_b);

  print(stdout, fmt, "Bond1", gpu::get_ebond());
}
}
