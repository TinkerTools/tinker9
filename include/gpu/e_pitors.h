#ifndef TINKER_GPU_E_PITORS_H_
#define TINKER_GPU_E_PITORS_H_

#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
extern int npitors;
extern int (*ipit)[6];
extern real* kpit;
extern real ptorunit;

extern real* ept;
extern real* vir_ept;

void epitors_data(rc_op op);

void epitors(int vers);
TINKER_NAMESPACE_END

#endif
