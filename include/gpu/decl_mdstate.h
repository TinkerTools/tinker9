#ifndef TINKER_GPU_DECL_MDSTATE_H_
#define TINKER_GPU_DECL_MDSTATE_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
enum {
  use_xyz = 0x001,   /// xyz
  use_vel = 0x002,   /// velocity
  use_accel = 0x004, /// acceleration
  use_mass = 0x008,  /// mass

  use_energy = 0x010, /// energy 16
  use_grad = 0x020,   /// gradient 32
  use_virial = 0x040, /// virial 64
  use_analyz = 0x080  /// analyze 128
};

template <int USE>
void sanity_check() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_a = USE & use_analyz;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  static_assert(do_v ? do_g : true, "");
  static_assert(do_a ? do_e : true, "");
}

extern int use_data;
extern int n;

extern real *x, *y, *z;
extern real *vx, *vy, *vz;
extern real *ax, *ay, *az;
extern real* mass;

void n_data();
void xyz_data(int op);
void vel_data(int op);
void accel_data(int op);
void mass_data(int op);
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
namespace gpu {
// clang-format off
const int v0 = use_energy;                         ///  16
const int v1 = use_energy + use_grad + use_virial; /// 112
const int v3 = use_energy + use_analyz;            /// 144

const int v4 = use_energy + use_grad;              ///  48
const int v5 = use_grad;                           ///  32
const int v6 = use_grad + use_virial;              ///  96
// clang-format on

extern real* esum;
extern real* vir;
extern real *gx, *gy, *gz;
const int _x = 0;
const int _y = 1;
const int _z = 2;
const int _xx = 0;
const int _yx = 1;
const int _zx = 2;
const int _xy = 3;
const int _yy = 4;
const int _zy = 5;
const int _xz = 6;
const int _yz = 7;
const int _zz = 8;

void zero_egv();
void egv_data(int op);

double get_energy(const real* e_gpu);
int get_count(const int* ecount_gpu);
void get_virial(double* v_out, const real* v_gpu);
}

/**
 * To solve the notorious trailing comma problem, see:
 * https://gcc.gnu.org/onlinedocs/cpp/Variadic-Macros.html
 */
#define TINKER_BONDED_DEF_1_(func, suffix, tmpl, vers, ...)                    \
  void func##suffix() { tmpl<vers, ##__VA_ARGS__>(); }
#define TINKER_BONDED_DEF(func, tmpl, ...)                                     \
  TINKER_BONDED_DEF_1_(func, 0, tmpl, gpu::v0, ##__VA_ARGS__)                  \
  TINKER_BONDED_DEF_1_(func, 1, tmpl, gpu::v1, ##__VA_ARGS__)                  \
  TINKER_BONDED_DEF_1_(func, 4, tmpl, gpu::v4, ##__VA_ARGS__)                  \
  TINKER_BONDED_DEF_1_(func, 5, tmpl, gpu::v5, ##__VA_ARGS__)                  \
  TINKER_BONDED_DEF_1_(func, 6, tmpl, gpu::v6, ##__VA_ARGS__)
#define TINKER_NONBONDED_DEF(func, tmpl, ...)                                  \
  TINKER_BONDED_DEF_1_(func, 0, tmpl, gpu::v0, ##__VA_ARGS__)                  \
  TINKER_BONDED_DEF_1_(func, 1, tmpl, gpu::v1, ##__VA_ARGS__)                  \
  TINKER_BONDED_DEF_1_(func, 3, tmpl, gpu::v3, ##__VA_ARGS__)                  \
  TINKER_BONDED_DEF_1_(func, 4, tmpl, gpu::v4, ##__VA_ARGS__)                  \
  TINKER_BONDED_DEF_1_(func, 5, tmpl, gpu::v5, ##__VA_ARGS__)                  \
  TINKER_BONDED_DEF_1_(func, 6, tmpl, gpu::v6, ##__VA_ARGS__)
TINKER_NAMESPACE_END

#endif
