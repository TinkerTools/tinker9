#ifndef TINKER_GPU_E_OPBEND_H_
#define TINKER_GPU_E_OPBEND_H_

#include "util_cxx.h"
#include "util_rc_man.h"

TINKER_NAMESPACE_BEGIN
typedef enum { opbend_w_d_c, opbend_allinger } eopbend_t;
extern eopbend_t opbtyp;

extern int nopbend;
extern int* iopb;
extern real* opbk;
extern real opbunit;
extern real copb, qopb, popb, sopb;

extern real* eopb;
extern real* vir_eopb;

void eopbend_data(rc_t rc);

void eopbend(int vers);
TINKER_NAMESPACE_END

#endif
