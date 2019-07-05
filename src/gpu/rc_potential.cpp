#include "gpu/decl_mdstate.h"
#include "gpu/decl_potent.h"
#include "gpu/e_potential.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
/// This function is used in void mdstate_data(rc_t).
void potential_data_(rc_t rc) {

  egv_data(rc);

  // bonded terms

  ebond_data(rc);
  eangle_data(rc);
  estrbnd_data(rc);
  eurey_data(rc);
  eopbend_data(rc);
  etors_data(rc);
  epitors_data(rc);
  etortor_data(rc);

  // non-bonded terms

  evdw_data(rc);

  // Must call elec_data() before any electrostatics routine.

  elec_data(rc);
  empole_data(rc);
  if (use_potent(polar_term))
    polargroup_data(rc);
  epolar_data(rc);
}

void energy_potential(int vers) {

  zero_egv(vers);

  // bonded terms

  if (use_potent(bond_term))
    ebond(vers);
  if (use_potent(angle_term))
    eangle(vers);
  if (use_potent(strbnd_term))
    estrbnd(vers);
  if (use_potent(urey_term))
    eurey(vers);
  if (use_potent(opbend_term))
    eopbend(vers);
  if (use_potent(torsion_term))
    etors(vers);
  if (use_potent(pitors_term))
    epitors(vers);
  if (use_potent(tortor_term))
    etortor(vers);

  // non-bonded terms

  if (use_potent(vdw_term))
    evdw(vers);

  gpu::elec_init(vers);
  if (use_potent(mpole_term))
    empole(vers);
  if (use_potent(polar_term))
    epolar(vers);
  gpu::torque(vers);

  sum_energies(vers);
}
}
TINKER_NAMESPACE_END
