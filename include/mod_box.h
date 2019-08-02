#ifndef TINKER_MOD_BOX_H_
#define TINKER_MOD_BOX_H_

#include "util_rc_man.h"

TINKER_NAMESPACE_BEGIN
/// @brief
/// periodic boundary conditions (pbc)
struct box_t {
  typedef enum {
    null = 0x000,  ///< no pbc
    ortho = 0x001, ///< orthogonal
    mono = 0x002,  ///< monoclinic
    tri = 0x004,   ///< triclinic
    oct = 0x008    ///< truncated octahedron
  } shape_t;

  /**
   * @brief
   * the matrix form of @c lvec and @c recip pbc box vectors
   *
   * @code
   * Triclinic
   * a.x,b.x,c.x = a.x b.x c.x
   * a.y,b.y,c.y =   0 b.y c.y
   * a.z,b.z,c.z =   0   0 c.z
   * @endcode
   *
   * @code
   * Monoclinic
   * alpha = gamma = 90 degrees
   * a.x,b.x,c.x = a.x   0 c.x
   * a.y,b.y,c.y =   0 b.y   0
   * a.z,b.z,c.z =   0   0 c.z
   * @endcode
   *
   * @code
   * Orthogonal
   * alpha = beta = gamma = 90 degrees
   * a.x,b.x,c.x = a.x   0   0
   * a.y,b.y,c.y =   0 b.y   0
   * a.z,b.z,c.z =   0   0 c.z
   * @endcode
   */
  /// @{
  real lvec[3][3];
  real recip[3][3];
  /// @}

  real volbox;   ///< volume of the pbc box
  shape_t shape; ///< shape of the pbc box
};

/// device pointer to the pbc box
TINKER_EXTERN box_t* box;
/// device pointer to the current pbc box of a trajectory
TINKER_EXTERN box_t* trajbox;

void box_data(rc_op);
void box_data_copyout(const box_t&);
TINKER_NAMESPACE_END

#endif
