#ifndef TINKER_GPU_E_PITORS_H_
#define TINKER_GPU_E_PITORS_H_

#include "util_cxx.h"
#include "util_rc_man.h"

TINKER_NAMESPACE_BEGIN
extern int npitors;
extern int (*ipit)[6];
extern real* kpit;
extern real ptorunit;

extern real* ept;
extern real* vir_ept;

void epitors_data(rc_t rc);

void epitors(int vers);
TINKER_NAMESPACE_END

#endif
