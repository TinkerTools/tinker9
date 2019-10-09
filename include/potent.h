#ifndef TINKER_POTENT_H_
#define TINKER_POTENT_H_

#include "macro.h"

TINKER_NAMESPACE_BEGIN
typedef enum
{
   bond_term,    // eb
   angle_term,   // ea
   strbnd_term,  // eba
   urey_term,    // eub
   angang_term,  // eaa
   opbend_term,  // eopb
   opdist_term,  // eopd
   imporp_term,  // eid
   imptors_term, // eit
   torsion_term, // et
   pitors_term,  // ept
   strtor_term,  // ebt
   angtor_term,  // eat
   tortor_term,  // ett

   vdw_term,    // ev
   repuls_term, // er
   disp_term,   // edsp
   charge_term, // ec
   chgdpl_term, // ecd
   dipole_term, // ed
   mpole_term,  // em
   polar_term,  // ep
   chgtrn_term, // ect
   rxnfld_term, // erxf
   solv_term,   // es
   metal_term,  // elf
   geom_term,   // eg
   extra_term   // ex
} potent_t;

int use_potent(potent_t term);
int count_bonded_term(potent_t term);
TINKER_NAMESPACE_END

#endif
