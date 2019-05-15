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

const int mpl_c = 0;
const int mpl_dx = 1;
const int mpl_dy = 2;
const int mpl_dz = 3;
const int mpl_qxx = 4;
const int mpl_qxy = 5;
const int mpl_qxz = 6;
const int mpl_qyx = mpl_qxy;
const int mpl_qyy = 7;
const int mpl_qyz = 8;
const int mpl_qzx = mpl_qxz;
const int mpl_qzy = mpl_qyz;
const int mpl_qzz = 9;
const int mpl_total = 10;
// (1), (2, 3, 4), (5, 6, 7, 9, 10, 13)
const int mpl_tinker[mpl_total] = {0, 1, 2, 3, 4, 5, 6, 8, 9, 12};
extern real (*pole)[mpl_total];
extern real (*rpole)[mpl_total];

extern real* em;
extern int* nem;
extern real* vir_em;
extern real* torque;

int use_empole();
void get_empole_type(int& typ, std::string& typ_str);
void e_mpole_data(int op);

void rotpole();
}
TINKER_NAMESPACE_END

#endif
