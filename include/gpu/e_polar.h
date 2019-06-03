#ifndef TINKER_GPU_E_POLAR_H_
#define TINKER_GPU_E_POLAR_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
extern real* polarity;
extern real* pdamp;
extern real* polarity_inv;

extern real* ep;
extern int* nep;
extern real* vir_ep;

int use_epolar();
void e_polar_data(int op);
}
TINKER_NAMESPACE_END

#endif
