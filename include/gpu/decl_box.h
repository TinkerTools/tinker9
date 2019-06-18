#ifndef TINKER_GPU_DECL_BOX_H_
#define TINKER_GPU_DECL_BOX_H_

#include "decl_real.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
enum {
  box_null = 0x000,  /// null
  box_ortho = 0x001, /// orthogonal
  box_mono = 0x002,  /// monoclinic
  box_tri = 0x004,   /// triclinic
  box_oct = 0x008    /// truncated octahedron
};

typedef struct box_def_st__ {
  real xbox, ybox, zbox;
  real alpha, beta, gamma;
  real lvec[3][3];
  real recip[3][3];
  real volbox;
  int shape;
} box_t;

extern box_t* box;

void box_data(int op);
}
TINKER_NAMESPACE_END

#endif
