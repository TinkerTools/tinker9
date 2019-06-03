#ifndef TINKER_GPU_E_POTENTIAL_H_
#define TINKER_GPU_E_POTENTIAL_H_

#include "e_bond.h"
#include "e_mpole.h"
#include "e_polar.h"
#include "e_vdw.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void potential_data(int op);
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_gradient1();
}

#endif
