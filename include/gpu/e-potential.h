#ifndef TINKER_GPU_E_POTENTIAL_H_
#define TINKER_GPU_E_POTENTIAL_H_

#include "e-bond.h"
#include "e-vdw.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void potential_data(int op);
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_gradient1();
}

#endif
