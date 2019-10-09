#ifndef TINKER_E_STRBND_H_
#define TINKER_E_STRBND_H_

#include "dev_array.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int nstrbnd;
TINKER_EXTERN device_pointer<int, 3> isb;
TINKER_EXTERN device_pointer<real, 2> sbk;
TINKER_EXTERN real stbnunit;

TINKER_EXTERN energy_buffer eba;
TINKER_EXTERN virial_buffer vir_eba;

void estrbnd_data(rc_op op);

void estrbnd(int vers);
TINKER_NAMESPACE_END

#endif
