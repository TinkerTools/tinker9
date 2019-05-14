#ifndef TINKER_GPU_DECL_BOX_H_
#define TINKER_GPU_DECL_BOX_H_

#include "decl-real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
const int box_null = 0x000;  /// null
const int box_ortho = 0x001; /// orthogonal
const int box_mono = 0x002;  /// monoclinic
const int box_tri = 0x004;   /// triclinic
const int box_oct = 0x008;   /// truncated octahedron

struct box_st {
  real xbox, ybox, zbox;
  real alpha, beta, gamma;
  real lvec[3][3];
  real recip[3][3];
  int shape;
};

extern box_st* box;

void box_data(int op);
}
TINKER_NAMESPACE_END

#endif
