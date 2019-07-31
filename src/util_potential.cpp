#include "util_potential.h"
#include "mod_md.h"
#include "mod_nblist.h"
#include "mod_polgrp.h"
#include "util_potent.h"

TINKER_NAMESPACE_BEGIN
void potential_data(rc_t rc) {
  if ((use_data & calc::vmask) == 0)
    return;

  rc_man egv42_{egv_data, rc};

  // bonded terms

  rc_man ebond42_{ebond_data, rc};
  rc_man eangle42_{eangle_data, rc};
  rc_man estrbnd42_{estrbnd_data, rc};
  rc_man eurey42_{eurey_data, rc};
  rc_man eopbend42_{eopbend_data, rc};
  rc_man etors42_{etors_data, rc};
  rc_man epitors42_{epitors_data, rc};
  rc_man etortor42_{etortor_data, rc};

  // non-bonded terms

  rc_man evdw42_{evdw_data, rc};

  // Must call elec_data() before any electrostatics routine.

  rc_man elec42_{elec_data, rc};

  rc_man empole42_{empole_data, rc};
  if (use_potent(polar_term))
    rc_man polargroup42_{polargroup_data, rc};
  rc_man epolar42_{epolar_data, rc};
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

  elec_init(vers);
  if (use_potent(mpole_term))
    empole(vers);
  if (use_potent(polar_term))
    epolar(vers);
  torque(vers);

  sum_energies(vers);

  // list update

  nblist_data(rc_evolve);
}
TINKER_NAMESPACE_END
