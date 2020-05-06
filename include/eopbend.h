#pragma once
#include "darray.h"
#include "energy_buffer.h"
#include "rc_man.h"

namespace tinker {
enum class eopbend_t
{
   w_d_c,
   allinger
};
TINKER_EXTERN eopbend_t opbtyp;

TINKER_EXTERN int nopbend;
TINKER_EXTERN pointer<int> iopb;
TINKER_EXTERN pointer<real> opbk;
TINKER_EXTERN real opbunit;
TINKER_EXTERN real copb, qopb, popb, sopb;

TINKER_EXTERN energy_buffer eopb;
TINKER_EXTERN virial_buffer vir_eopb;
TINKER_EXTERN grad_prec *deopbx, *deopby, *deopbz;
TINKER_EXTERN energy_prec energy_eopb;

void eopbend_data(rc_op op);

void eopbend(int vers);
void eopbend_acc(int);
}
