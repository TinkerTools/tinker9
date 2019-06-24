#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_nblist.h"
#include "gpu/decl_pme.h"
#include "gpu/decl_potent.h"
#include "gpu/e_potential.h"
#include "util/format_print.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void potential_data(int op) {
  ebond_data(op);
  eangle_data(op);
  estrbnd_data(op);

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

  if (use_potent(bond_term)) {
    ebond(vers);
  }

  if (use_evdw()) {
    evdw(vers);
  }

  print(stdout, title, "Kcal/mole", "Interactions");

  if (use_potent(bond_term))
    print(stdout, fmt, "Bond Stretching", get_energy(eb),
          count_bonded_term(bond_term), bndtyp_str);

  if (use_evdw())
    print(stdout, fmt, "Van der Waals", get_energy(ev), get_count(nev),
          vdwtyp_str);
}
}
TINKER_NAMESPACE_END
