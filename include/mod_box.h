#ifndef TINKER_MOD_BOX_H_
#define TINKER_MOD_BOX_H_

#include "util_rc_man.h"

TINKER_NAMESPACE_BEGIN
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
struct box_t {
  typedef enum {
    null = 0x000,  /// null
    ortho = 0x001, /// orthogonal
    mono = 0x002,  /// monoclinic
    tri = 0x004,   /// triclinic
    oct = 0x008    /// truncated octahedron
  } shape_t;

  real lvec[3][3];
  real recip[3][3];
  real volbox;
  shape_t shape;
};

TINKER_EXTERN box_t* box;
TINKER_EXTERN box_t* trajbox;

void box_data(rc_t);
void box_data_copyout(const box_t&);
TINKER_NAMESPACE_END

#endif
