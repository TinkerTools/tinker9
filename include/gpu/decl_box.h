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

/**
 * matrix form of lvec/recip periodic box vectors
 *
 * triclinic
 * a.x,b.x,c.x = a.x b.x c.x
 * a.y,b.y,c.y =   0 b.y c.y
 * a.z,b.z,c.z =   0   0 c.z
 *
 * monoclinic alpha = gamma = 90
 * a.x,b.x,c.x = a.x   0 c.x
 * a.y,b.y,c.y =   0 b.y   0
 * a.z,b.z,c.z =   0   0 c.z
 *
 * orthogonal alpha = beta = gamma = 90
 * a.x,b.x,c.x = a.x   0   0
 * a.y,b.y,c.y =   0 b.y   0
 * a.z,b.z,c.z =   0   0 c.z
 *
 * cartesian_column_vector dot recip = fractional_column_vector
 * fractional_column_vector dot lvec = cartesian_column_vector
 */
typedef struct box_def_st__ {
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
