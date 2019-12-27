#include "energy.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"

TINKER_NAMESPACE_BEGIN
void potential_data(rc_op op)
{
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

   // misc. terms

   rc_man egeom42_{egeom_data, op};

   // non-bonded terms

   rc_man evdw42_{evdw_data, op};

   // Must call elec_data() before any electrostatics routine.

   rc_man elec42_{elec_data, op};

   rc_man empole42_{empole_data, op};
   rc_man epolar42_{epolar_data, op};
#if TINKER_CUDART
   // Must follow empole_data() and epolar_data().
   rc_man emplar42_{emplar_data, op};
#endif
}

void energy_potential(int vers)
{

   zero_egv(vers);

   wait_queue();

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

   // misc. terms

   if (use_potent(geom_term))
      egeom(vers);

   // non-bonded terms

   if (use_potent(vdw_term))
      evdw(vers);

   elec_init(vers);
#if TINKER_CUDART
   if (use_potent(mpole_term) && use_potent(polar_term) &&
       !(vers & calc::analyz) && mlist_version() == NBList::spatial) {
      emplar_cu(vers);
      goto skip_mpole_polar;
   }
#endif
   if (use_potent(mpole_term))
      empole(vers);
   if (use_potent(polar_term))
      epolar(vers);
skip_mpole_polar:
   torque(vers);
   wait_queue();

   sum_energies(vers);
}
TINKER_NAMESPACE_END
