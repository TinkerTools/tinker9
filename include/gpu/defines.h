#ifndef TINKER_GPU_DEFINES_H_
#define TINKER_GPU_DEFINES_H_

#include "mathfunc.h"
#include <string>

TINKER_NAMESPACE_BEGIN
namespace gpu {
#ifdef TINKER_GPU_DOUBLE
typedef double real;
#endif

#ifdef TINKER_GPU_SINGLE
typedef float real;
#endif

const int op_destroy = 0;
const int op_create = 1;

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

const int box_null  = 0x000; /// null
const int box_ortho = 0x001; /// orthogonal
const int box_mono  = 0x002; /// monoclinic
const int box_tri   = 0x004; /// triclinic
const int box_oct   = 0x008; /// truncated octahedron
// clang-format on

namespace gpu {
struct box_st {
  real xbox, ybox, zbox;
  real alpha, beta, gamma;
  real lvec[3][3];
  real recip[3][3];
  int shape;
};

struct couple_st {
  static const int maxn12 = 8;
  static const int maxn13 = 32;
  static const int maxn14 = 96;  // 128;
  static const int maxn15 = 224; // 256;

  int *n12, *n13, *n14, *n15;
  int (*i12)[maxn12];
  int (*i13)[maxn13];
  int (*i14)[maxn14];
  int (*i15)[maxn15];
};

struct nblist_st {
  int* nlst;
  int* lst;
  real *xold, *yold, *zold;
  const real* x;
  const real* y;
  const real* z;
  int maxnlst;
  real cutoff, buffer;
};
}
TINKER_NAMESPACE_END

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
#endif
