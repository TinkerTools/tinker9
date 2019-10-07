#ifndef TINKER_E_MPOLE_H_
#define TINKER_E_MPOLE_H_

#include "elec.h"
#include "energy_buffer.h"

TINKER_NAMESPACE_BEGIN
/// \defgroup mpole Multipole Electrostatic Energy
/// \ingroup gvar

TINKER_EXTERN elec_t empole_electyp;

TINKER_EXTERN real m2scale, m3scale, m4scale, m5scale;

TINKER_EXTERN int nmexclude_;
TINKER_EXTERN device_pointer<int, 2> mexclude_;
TINKER_EXTERN device_pointer<real> mexclude_scale_;

TINKER_EXTERN NonbondedEnergy em_handle;

void empole_data (rc_op op);

void empole_coulomb (int vers);
void empole_ewald (int vers);
void empole (int vers);
TINKER_NAMESPACE_END

#endif
