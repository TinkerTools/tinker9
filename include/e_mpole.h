#ifndef TINKER_E_MPOLE_H_
#define TINKER_E_MPOLE_H_

#include "elec.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN elec_t empole_electyp;

TINKER_EXTERN real m2scale, m3scale, m4scale, m5scale;

TINKER_EXTERN real* em;
TINKER_EXTERN int* nem;
TINKER_EXTERN real* vir_em;

void empole_data(rc_op op);

void empole_coulomb(int vers);
void empole_ewald(int vers);
void empole(int vers);
TINKER_NAMESPACE_END

#endif
