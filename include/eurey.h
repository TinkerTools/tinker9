#pragma once
#include "darray.h"
#include "energy_buffer.h"
#include "rc_man.h"

namespace tinker {
// module urey
TINKER_EXTERN int nurey;
TINKER_EXTERN pointer<int, 3> iury;
TINKER_EXTERN pointer<real> uk, ul;

// module urypot
TINKER_EXTERN real cury, qury, ureyunit;

void eurey_data(rc_op op);

void eurey(int vers);
void eurey_acc(int);
}
