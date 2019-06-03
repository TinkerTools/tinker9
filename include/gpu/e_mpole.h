#ifndef TINKER_GPU_E_MPOLE_H_
#define TINKER_GPU_E_MPOLE_H_

#include "decl_elec.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
extern int electyp;
extern std::string electyp_str;

extern double mpole_switch_cut, mpole_switch_off;

extern real* em;
extern int* nem;
extern real* vir_em;

int use_empole();
void get_empole_type(int& typ, std::string& typ_str);
void e_mpole_data(int op);

/**
 * @brief
 * "grid_mpole" places the fractional atomic multipoles onto
 * the particle mesh Ewald grid.
 */
void grid_mpole(real (*gpu_fmp)[10]);

/**
 * @brief
 * "fphi_mpole" extracts the permanent multipole potential from
 * the particle mesh Ewald grid.
 */
void fphi_mpole(real (*gpu_fphi)[20]);

/**
 * @brief
 * make the scalar summation over reciprocal lattice
 */
void pme_conv0(int pme_unit);                 // without virial
void pme_conv1(int pme_unit, real* gpu_vir9); // with virial
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_empole_coulomb0();
void tinker_gpu_empole_coulomb1();
void tinker_gpu_empole_coulomb3();
void tinker_gpu_empole_coulomb4();
void tinker_gpu_empole_coulomb5();
void tinker_gpu_empole_coulomb6();

void tinker_gpu_emreal0();
void tinker_gpu_emreal1();
void tinker_gpu_emreal3();
void tinker_gpu_emreal4();
void tinker_gpu_emreal5();
void tinker_gpu_emreal6();

void tinker_gpu_emrecip0();
void tinker_gpu_emrecip1();
void tinker_gpu_emrecip3();
void tinker_gpu_emrecip4();
void tinker_gpu_emrecip5();
void tinker_gpu_emrecip6();

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
