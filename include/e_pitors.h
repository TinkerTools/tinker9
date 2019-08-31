#ifndef TINKER_E_PITORS_H_
#define TINKER_E_PITORS_H_

#include "device_vector.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int npitors;
TINKER_EXTERN DeviceVector<int, 6> ipit_vec;
TINKER_EXTERN DeviceVector<real> kpit_vec;
TINKER_EXTERN real ptorunit;

TINKER_EXTERN BondedEnergy ept_handle;

void epitors_data(rc_op op);

void epitors(int vers);
TINKER_NAMESPACE_END

#endif
