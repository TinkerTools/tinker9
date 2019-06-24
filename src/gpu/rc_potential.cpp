#include "gpu/decl_mdstate.h"
#include "gpu/decl_nblist.h"
#include "gpu/decl_pme.h"
#include "gpu/decl_potent.h"
#include "gpu/e_potential.h"
#include "util/format_print.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void potential_data(rc_t rc) {
  ebond_data(rc);
  eangle_data(rc);
  estrbnd_data(rc);

  evdw_data(rc);

  // Must call elec_data() before any electrostatics routine.

  elec_data(rc);
  empole_data(rc);
  if (use_epolar())
    polargroup_data(rc);
  epolar_data(rc);
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
