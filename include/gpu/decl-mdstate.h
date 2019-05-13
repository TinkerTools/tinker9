#ifndef TINKER_GPU_DECL_MDSTATE_H_
#define TINKER_GPU_DECL_MDSTATE_H_

#include "decl-box.h"
#include "decl-couple.h"
#include "decl-nblist.h"

TINKER_NAMESPACE_BEGIN
// clang-format off
const int use_xyz    = 0x001; /// xyz
const int use_vel    = 0x002; /// velocity
const int use_accel  = 0x004; /// acceleration
const int use_mass   = 0x008; /// mass
const int use_energy = 0x010; /// energy 16
const int use_grad   = 0x020; /// gradient 32
const int use_virial = 0x040; /// virial 64
const int use_analyz = 0x080; /// analyze 128

namespace gpu {
const int v0 = use_energy;                         ///  16
const int v1 = use_energy + use_grad + use_virial; /// 112
const int v3 = use_energy + use_analyz;            /// 144

const int v4 = use_energy + use_grad;              ///  48
const int v5 = use_grad;                           ///  32
const int v6 = use_grad + use_virial;              ///  96
}
// clang-format on

/**
 * To solve the notorious trailing comma problem, see:
 * https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
 */
#define TINKER_BONDED_GEN_1_(func, suffix, tmpl, vers, ...)                    \
  void func##suffix() { tmpl<vers, ##__VA_ARGS__>(); }
#define TINKER_BONDED_GEN(func, tmpl, ...)                                     \
  TINKER_BONDED_GEN_1_(func, 0, tmpl, gpu::v0, ##__VA_ARGS__)                  \
  TINKER_BONDED_GEN_1_(func, 1, tmpl, gpu::v1, ##__VA_ARGS__)                  \
  TINKER_BONDED_GEN_1_(func, 4, tmpl, gpu::v4, ##__VA_ARGS__)                  \
  TINKER_BONDED_GEN_1_(func, 5, tmpl, gpu::v5, ##__VA_ARGS__)                  \
  TINKER_BONDED_GEN_1_(func, 6, tmpl, gpu::v6, ##__VA_ARGS__)
#define TINKER_NONBONDED_GEN(func, tmpl, ...)                                  \
  TINKER_BONDED_GEN_1_(func, 0, tmpl, gpu::v0, ##__VA_ARGS__)                  \
  TINKER_BONDED_GEN_1_(func, 1, tmpl, gpu::v1, ##__VA_ARGS__)                  \
  TINKER_BONDED_GEN_1_(func, 3, tmpl, gpu::v3, ##__VA_ARGS__)                  \
  TINKER_BONDED_GEN_1_(func, 4, tmpl, gpu::v4, ##__VA_ARGS__)                  \
  TINKER_BONDED_GEN_1_(func, 5, tmpl, gpu::v5, ##__VA_ARGS__)                  \
  TINKER_BONDED_GEN_1_(func, 6, tmpl, gpu::v6, ##__VA_ARGS__)

namespace gpu {
extern int use_data;
extern int n;

extern real* esum;
extern real* vir;
extern real* mass;

extern real *x, *y, *z;
extern real *vx, *vy, *vz;
extern real *ax, *ay, *az;
extern real *gx, *gy, *gz;

const int _xx = 0;
const int _yx = 1;
const int _zx = 2;
const int _xy = 3;
const int _yy = 4;
const int _zy = 5;
const int _xz = 6;
const int _yz = 7;
const int _zz = 8;
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_zero_vag();
}

#endif
