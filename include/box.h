#ifndef TINKER_BOX_H_
#define TINKER_BOX_H_

#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
/// @brief
/// periodic boundary conditions (pbc)
struct Box {
  typedef enum {
    null = 0x000,  ///< no pbc
    ortho = 0x001, ///< orthogonal
    mono = 0x002,  ///< monoclinic
    tri = 0x004,   ///< triclinic
    oct = 0x008    ///< truncated octahedron
  } Shape;

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
   * Monoclinic: alpha = gamma = 90 degrees
   * a.x,b.x,c.x = a.x   0 c.x
   * a.y,b.y,c.y =   0 b.y   0
   * a.z,b.z,c.z =   0   0 c.z
   * @endcode
   *
   * @code
   * Orthogonal: alpha = beta = gamma = 90 degrees
   * a.x,b.x,c.x = a.x   0   0
   * a.y,b.y,c.y =   0 b.y   0
   * a.z,b.z,c.z =   0   0 c.z
   * @endcode
   */
  /// @{
  real lvec[3][3];
  real recip[3][3];
  /// @}

  real volbox; ///< volume of the pbc box
  Shape shape; ///< shape of the pbc box
};

/// device pointer to the pbc box
TINKER_EXTERN Box* box;
/// device pointer to the current pbc box of a trajectory
TINKER_EXTERN Box* trajbox;

void box_data(rc_op);
void copyout_box_data(const Box*);
TINKER_NAMESPACE_END

#endif
