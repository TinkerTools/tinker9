#ifndef TINKER_E_UREY_H_
#define TINKER_E_UREY_H_

#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
// module urey
TINKER_EXTERN int nurey;
TINKER_EXTERN int (*iury)[3];
TINKER_EXTERN real *uk, *ul;

// module urypot
TINKER_EXTERN real cury, qury, ureyunit;

TINKER_EXTERN BondedEnergy eub_handle;

void eurey_data(rc_op op);

void eurey(int vers);
TINKER_NAMESPACE_END

#endif
