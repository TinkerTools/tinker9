#ifndef TINKER_E_PITORS_H_
#define TINKER_E_PITORS_H_

#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int npitors;
TINKER_EXTERN int (*ipit)[6];
TINKER_EXTERN real* kpit;
TINKER_EXTERN real ptorunit;

TINKER_EXTERN real* ept;
TINKER_EXTERN real* vir_ept;

void epitors_data(rc_op op);

void epitors(int vers);
TINKER_NAMESPACE_END

#endif
