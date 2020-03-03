#pragma once
#include "darray.h"
#include "energy_buffer.h"
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
enum class ebond_t
{
   harmonic,
   morse
};
TINKER_EXTERN ebond_t bndtyp;


TINKER_EXTERN real cbnd, qbnd, bndunit;
TINKER_EXTERN int nbond;
TINKER_EXTERN pointer<int, 2> ibnd;
TINKER_EXTERN pointer<real> bl, bk;


TINKER_EXTERN energy_buffer eb;
TINKER_EXTERN virial_buffer vir_eb;


void ebond_data(rc_op op);
void ebond(int vers);
TINKER_NAMESPACE_END
