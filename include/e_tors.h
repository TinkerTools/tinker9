#ifndef TINKER_E_TORS_H_
#define TINKER_E_TORS_H_

#include "device_vector.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int ntors;
TINKER_EXTERN DeviceVector<int, 4> itors_vec;
TINKER_EXTERN DeviceVector<real, 4> tors1_vec, tors2_vec, tors3_vec, tors4_vec,
    tors5_vec, tors6_vec;
TINKER_EXTERN real torsunit;

TINKER_EXTERN BondedEnergy et_handle;

void etors_data(rc_op op);

void etors(int vers);
TINKER_NAMESPACE_END

#endif
