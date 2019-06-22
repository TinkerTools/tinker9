#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_nblist.h"
#include "gpu/decl_pme.h"
#include "gpu/e_potential.h"
#include "util/format_print.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void potential_data(int op) {
  ebond_data(op);

  evdw_data(op);

  // Must call elec_data() before any electrostatics routine.

  elec_data(op);
  empole_data(op);
  if (use_epolar())
    polargroup_data(op);
  epolar_data(op);
}

void gradient(int vers) {
  const char* title = " Energy Component Breakdown :{:>20s}{:>20s}\n\n";
  const char* fmt = " {:28s}{:>20.4f}{:>17d}        {}\n";

  if (gpu::use_ebond()) {
    gpu::ebond(vers);
  }

  if (gpu::use_evdw()) {
    gpu::evdw(vers);
  }

  print(stdout, title, "Kcal/mole", "Interactions");

  if (gpu::use_ebond())
    print(stdout, fmt, "Bond Stretching", gpu::get_energy(gpu::eb),
          gpu::count_ebond(), gpu::bndtyp_str);

  if (gpu::use_evdw())
    print(stdout, fmt, "Van der Waals", gpu::get_energy(gpu::ev),
          gpu::get_count(gpu::nev), gpu::vdwtyp_str);
}
}
TINKER_NAMESPACE_END
