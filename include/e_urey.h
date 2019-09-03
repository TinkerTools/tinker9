#ifndef TINKER_E_UREY_H_
#define TINKER_E_UREY_H_

#include "dev_array.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
// module urey
TINKER_EXTERN int nurey;
TINKER_EXTERN device_array::ptr<int, 3>::type iury;
TINKER_EXTERN device_array::ptr<real>::type uk, ul;

// module urypot
TINKER_EXTERN real cury, qury, ureyunit;

TINKER_EXTERN BondedEnergy eub_handle;

void eurey_data(rc_op op);

void eurey(int vers);
TINKER_NAMESPACE_END

#endif
