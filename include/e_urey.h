#ifndef TINKER_E_UREY_H_
#define TINKER_E_UREY_H_

#include "device_vector.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
// module urey
TINKER_EXTERN int nurey;
TINKER_EXTERN DeviceVector<int, 3> iury_vec;
TINKER_EXTERN DeviceVector<real> uk_vec, ul_vec;

// module urypot
TINKER_EXTERN real cury, qury, ureyunit;

TINKER_EXTERN BondedEnergy eub_handle;

void eurey_data(rc_op op);

void eurey(int vers);
TINKER_NAMESPACE_END

#endif
