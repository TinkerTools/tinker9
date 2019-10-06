#ifndef TINKER_E_PITORS_H_
#define TINKER_E_PITORS_H_

#include "dev_array.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int npitors;
TINKER_EXTERN device_pointer<int, 6> ipit;
TINKER_EXTERN device_pointer<real> kpit;
TINKER_EXTERN real ptorunit;

TINKER_EXTERN BondedEnergy ept_handle;

void epitors_data (rc_op op);

void epitors (int vers);
TINKER_NAMESPACE_END

#endif
