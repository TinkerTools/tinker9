#include "gpu/e-potential.h"
#include "gpu/decl-dataop.h"
#include "gpu/decl-nblist.h"
#include "util/format.print.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
double get_energy(const real* e_gpu) {
  double e_out;
  copyout_data_1(&e_out, e_gpu, 1);
  return e_out;
}

void get_virial(double* v_out, const real* v_gpu) {
  copyout_data_1(v_out, v_gpu, 9);
}

void potential_data(int op) {
  e_bond_data(op);
  e_vdw_data(op);
}
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_gradient1() {
  m_tinker_using_namespace;
  const char* title = " Energy Component Breakdown :{:>20s}{:>20s}\n\n";
  const char* fmt = " {:28s}{:>20.4f}{:>17d}        {}\n";

  if (gpu::use_ebond()) {
    tinker_gpu_ebond_harmonic0();
    tinker_gpu_ebond_harmonic4();
    tinker_gpu_ebond_harmonic5();
    tinker_gpu_ebond_harmonic6();
    tinker_gpu_ebond_harmonic1();
  }

  if (gpu::use_evdw()) {
    /*
    tinker_gpu_evdw0();
    tinker_gpu_evdw3();
    tinker_gpu_evdw4();
    tinker_gpu_evdw5();
    tinker_gpu_evdw6();
    tinker_gpu_evdw1();
    // */
    tinker_gpu_evdw3();
  }

  tinker_gpu_vlist_update();

  print(stdout, title, "Kcal/mole", "Interactions");

  if (gpu::use_ebond())
    print(stdout, fmt, "Bond Stretching", gpu::get_energy(gpu::eb),
          gpu::count_ebond(), gpu::bndtyp_str);

  if (gpu::use_evdw())
    print(stdout, fmt, "Van der Waals", gpu::get_energy(gpu::ev),
          gpu::count_evdw(), gpu::vdwtyp_str);
}
}
