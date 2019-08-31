#ifndef TINKER_E_BOND_H_
#define TINKER_E_BOND_H_

#include "device_vector.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
enum class ebond_t { harmonic, morse };
TINKER_EXTERN ebond_t bndtyp;

TINKER_EXTERN real cbnd, qbnd, bndunit;
TINKER_EXTERN int nbond;
TINKER_EXTERN DeviceVector<int, 2> ibnd_vec;
TINKER_EXTERN DeviceVector<real> bl_vec, bk_vec;

TINKER_EXTERN BondedEnergy eb_handle;

void ebond_data(rc_op op);

void ebond(int vers);
TINKER_NAMESPACE_END

#endif
