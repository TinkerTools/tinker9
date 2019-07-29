#ifndef TINKER_GPU_E_UREY_H_
#define TINKER_GPU_E_UREY_H_

#include "util_cxx.h"
#include "util_rc_man.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
// module urey
extern int nurey;
extern int (*iury)[3];
extern real *uk, *ul;

// module urypot
extern real cury, qury, ureyunit;

extern real* eub;
extern real* vir_eub;

void eurey_data(rc_t rc);

void eurey(int vers);
}
TINKER_NAMESPACE_END

#endif
