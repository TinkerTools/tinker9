#ifndef TINKER_GPU_E_BOND_H_
#define TINKER_GPU_E_BOND_H_

#include "defines.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
const int ebond_harmonic = 0x001;
const int ebond_morse = 0x002;
extern int bndtyp;
extern std::string bndtyp_str;

extern real cbnd, qbnd, bndunit;
extern int nbond;
extern int (*ibnd)[2];
extern real *bl, *bk;

extern real* eb;
int use_ebond();
real get_ebond();
int count_ebond();
void e_bond_data(int op);
}
TINKER_NAMESPACE_END

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

#endif
