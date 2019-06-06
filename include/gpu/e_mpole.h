#ifndef TINKER_GPU_E_MPOLE_H_
#define TINKER_GPU_E_MPOLE_H_

#include "decl_elec.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
extern int empole_electyp;
extern std::string empole_electyp_str;

extern double mpole_switch_cut, mpole_switch_off;

extern real* em;
extern int* nem;
extern real* vir_em;

int use_empole();
void get_empole_type(int& typ, std::string& typ_str);
void e_mpole_data(int op);
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_empole_coulomb0();
void tinker_gpu_empole_coulomb1();
void tinker_gpu_empole_coulomb3();
void tinker_gpu_empole_coulomb4();
void tinker_gpu_empole_coulomb5();
void tinker_gpu_empole_coulomb6();

void tinker_gpu_empole_ewald0();
void tinker_gpu_empole_ewald1();
void tinker_gpu_empole_ewald3();
void tinker_gpu_empole_ewald4();
void tinker_gpu_empole_ewald5();
void tinker_gpu_empole_ewald6();

void tinker_gpu_empole0();
void tinker_gpu_empole1();
void tinker_gpu_empole3();
void tinker_gpu_empole4();
void tinker_gpu_empole5();
void tinker_gpu_empole6();
}

#endif
