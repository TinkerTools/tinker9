#ifndef TINKER_E_STRBND_H_
#define TINKER_E_STRBND_H_

#include "device_vector.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int nstrbnd;
TINKER_EXTERN DeviceVector<int, 3> isb_vec;
TINKER_EXTERN DeviceVector<real, 2> sbk_vec;
TINKER_EXTERN real stbnunit;

TINKER_EXTERN BondedEnergy eba_handle;

void estrbnd_data(rc_op op);

void estrbnd(int vers);
TINKER_NAMESPACE_END

#endif
