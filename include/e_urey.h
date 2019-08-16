#ifndef TINKER_E_UREY_H_
#define TINKER_E_UREY_H_

#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
// module urey
TINKER_EXTERN int nurey;
TINKER_EXTERN int (*iury)[3];
TINKER_EXTERN real *uk, *ul;

// module urypot
TINKER_EXTERN real cury, qury, ureyunit;

TINKER_EXTERN real* eub;
TINKER_EXTERN real* vir_eub;

void eurey_data(rc_op op);

void eurey(int vers);
TINKER_NAMESPACE_END

#endif
