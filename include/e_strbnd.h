#ifndef TINKER_E_STRBND_H_
#define TINKER_E_STRBND_H_

#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int nstrbnd;
TINKER_EXTERN int (*isb)[3];
TINKER_EXTERN real (*sbk)[2];
TINKER_EXTERN real stbnunit;

TINKER_EXTERN BondedEnergy eba_handle;

void estrbnd_data(rc_op op);

void estrbnd(int vers);
TINKER_NAMESPACE_END

#endif
