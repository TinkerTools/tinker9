#pragma once
#include "darray.h"
#include "energy_buffer.h"
#include "rc_man.h"

namespace tinker {
TINKER_EXTERN int npitors;
TINKER_EXTERN pointer<int, 6> ipit;
TINKER_EXTERN pointer<real> kpit;
TINKER_EXTERN real ptorunit;

void epitors_data(rc_op op);

void epitors(int vers);
void epitors_acc(int);
}
