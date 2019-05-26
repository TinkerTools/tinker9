#ifndef TINKER_GPU_E_MPOLE_H_
#define TINKER_GPU_E_MPOLE_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
const int elec_coulomb = 1;
const int elec_ewald = 2;
extern int electyp;
extern std::string electyp_str;

const int pole_none = 0;
const int pole_z_only = 1;
const int pole_z_then_x = 2;
const int pole_bisector = 3;
const int pole_z_bisect = 4;
const int pole_3_fold = 5;
struct local_frame_def_st {
  int zaxis, xaxis, yaxis;
  int polaxe;
};
extern local_frame_def_st* zaxis;

// PME: 0, x, y, z, xx, yy, zz, xy, xz, yz
enum {
  mpl_pme_0 = 0,
  mpl_pme_x = 1,
  mpl_pme_y = 2,
  mpl_pme_z = 3,
  mpl_pme_xx = 4,
  mpl_pme_yy = 5,
  mpl_pme_zz = 6,
  mpl_pme_xy = 7,
  mpl_pme_xz = 8,
  mpl_pme_yz = 9,
  mpl_pme_yx = mpl_pme_xy,
  mpl_pme_zx = mpl_pme_xz,
  mpl_pme_zy = mpl_pme_yz,
};
const int mpl_total = 10;
extern real (*pole)[mpl_total];
extern real (*rpole)[mpl_total];

extern double mpole_switch_cut, mpole_switch_off;

extern real* em;
extern int* nem;
extern real* vir_em;
extern real *trqx, *trqy, *trqz;

int use_empole();
void get_empole_type(int& typ, std::string& typ_str);
void e_mpole_data(int op);

void chkpole();
void rotpole();
void zero_torque();
void torque0(); // gradient only
void torque1(); // gradient and virial

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
void fphi_mpole();
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
