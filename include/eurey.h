#pragma once
#include "darray.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
// module urey
TINKER_EXTERN int nurey;
TINKER_EXTERN pointer<int, 3> iury;
TINKER_EXTERN pointer<real> uk, ul;

// module urypot
TINKER_EXTERN real cury, qury, ureyunit;

TINKER_EXTERN energy_buffer eub;
TINKER_EXTERN virial_buffer vir_eub;
TINKER_EXTERN grad_prec *deubx, *deuby, *deubz;
TINKER_EXTERN energy_prec energy_eub;

void eurey_data(rc_op op);

void eurey(int vers);
void eurey_acc(int);
TINKER_NAMESPACE_END
