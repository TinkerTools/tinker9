#ifndef TINKER_GPU_DEFINES_H_
#define TINKER_GPU_DEFINES_H_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
#ifdef TINKER_GPU_DOUBLE
typedef double real;
#endif

#ifdef TINKER_GPU_SINGLE
typedef float real;
#endif
}

// clang-format off
const int use_xyz    = 0x001; /// xyz
const int use_vel    = 0x002; /// velocity
const int use_accel  = 0x004; /// acceleration
const int use_mass   = 0x008; /// mass
const int use_energy = 0x010; /// energy
const int use_grad   = 0x020; /// gradient
const int use_virial = 0x040; /// virial

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
}
TINKER_NAMESPACE_END

#endif
