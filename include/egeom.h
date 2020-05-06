#pragma once
#include "darray.h"
#include "energy_buffer.h"
#include "rc_man.h"


namespace tinker {
// group distance restraints
/// \ingroup geom
/// \brief Number of group distance restraints to be applied.
TINKER_EXTERN int ngfix;
/// \ingroup geom
/// \brief Group numbers defining each group distance restraint.
TINKER_EXTERN pointer<int, 2> igfix;
/// \ingroup geom
/// \brief Force constant and target range for each group distance.
TINKER_EXTERN pointer<real, 3> gfix;


TINKER_EXTERN energy_buffer eg;
TINKER_EXTERN virial_buffer vir_eg;
TINKER_EXTERN grad_prec *degx, *degy, *degz;
TINKER_EXTERN energy_prec energy_eg;


void egeom_data(rc_op);

void egeom(int vers);
void egeom_acc(int);
}
