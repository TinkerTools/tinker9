#include "energy.h"
#include "md.h"
#include "nblist.h"
#include "polgrp.h"
#include "potent.h"

TINKER_NAMESPACE_BEGIN
void potential_data(rc_op op) {
  if ((rc_flag & calc::vmask) == 0)
    return;

  rc_man egv42_{egv_data, op};

  // bonded terms

  rc_man ebond42_{ebond_data, op};
  rc_man eangle42_{eangle_data, op};
  rc_man estrbnd42_{estrbnd_data, op};
  rc_man eurey42_{eurey_data, op};
  rc_man eopbend42_{eopbend_data, op};
  rc_man etors42_{etors_data, op};
  rc_man epitors42_{epitors_data, op};
  rc_man etortor42_{etortor_data, op};

  // non-bonded terms

  rc_man evdw42_{evdw_data, op};

  // Must call elec_data() before any electrostatics routine.

  rc_man elec42_{elec_data, op};

  rc_man empole42_{empole_data, op};
  if (use_potent(polar_term))
    rc_man polargroup42_{polargroup_data, op};
  rc_man epolar42_{epolar_data, op};
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

  nblist_data(rc_man::evolve);
}
TINKER_NAMESPACE_END
