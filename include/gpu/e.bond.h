#ifndef TINKER_GPU_E_BOND_H_
#define TINKER_GPU_E_BOND_H_

#include "defines.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void e_bond_data(int op);

extern "C" {
void tinker_gpu_ebond_harmonic0();
void tinker_gpu_ebond_harmonic1();
void tinker_gpu_ebond_harmonic4();
void tinker_gpu_ebond_harmonic5();
void tinker_gpu_ebond_harmonic6();

void tinker_gpu_ebond_morse0();
void tinker_gpu_ebond_morse1();
void tinker_gpu_ebond_morse4();
void tinker_gpu_ebond_morse5();
void tinker_gpu_ebond_morse6();
}
}
TINKER_NAMESPACE_END

#endif
