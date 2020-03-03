#pragma once
#include "dev_array.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int npitors;
TINKER_EXTERN pointer<int, 6> ipit;
TINKER_EXTERN pointer<real> kpit;
TINKER_EXTERN real ptorunit;

TINKER_EXTERN energy_buffer ept;
TINKER_EXTERN virial_buffer vir_ept;

void epitors_data(rc_op op);

void epitors(int vers);
TINKER_NAMESPACE_END
