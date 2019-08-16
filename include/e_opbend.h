#ifndef TINKER_E_OPBEND_H_
#define TINKER_E_OPBEND_H_

#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
enum class eopbend_t { w_d_c, allinger };
TINKER_EXTERN eopbend_t opbtyp;

TINKER_EXTERN int nopbend;
TINKER_EXTERN int* iopb;
TINKER_EXTERN real* opbk;
TINKER_EXTERN real opbunit;
TINKER_EXTERN real copb, qopb, popb, sopb;

TINKER_EXTERN real* eopb;
TINKER_EXTERN real* vir_eopb;

void eopbend_data(rc_op op);

void eopbend(int vers);
TINKER_NAMESPACE_END

#endif
