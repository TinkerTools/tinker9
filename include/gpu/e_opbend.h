#ifndef TINKER_GPU_E_OPBEND_H_
#define TINKER_GPU_E_OPBEND_H_

#include "decl_basic.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
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

void eopbend_acc_impl__(int vers);
void eopbend(int vers);
}
TINKER_NAMESPACE_END

#endif
