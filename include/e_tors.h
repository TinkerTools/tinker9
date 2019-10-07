#ifndef TINKER_E_TORS_H_
#define TINKER_E_TORS_H_

#include "dev_array.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int ntors;
TINKER_EXTERN device_pointer<int, 4> itors;
TINKER_EXTERN device_pointer<real, 4> tors1, tors2, tors3, tors4, tors5, tors6;
TINKER_EXTERN real torsunit;

TINKER_EXTERN BondedEnergy et_handle;

void etors_data (rc_op op);

void etors (int vers);
TINKER_NAMESPACE_END

#endif
