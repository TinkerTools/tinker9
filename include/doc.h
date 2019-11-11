/**
 * \defgroup ff      Force Field Terms
 * \defgroup io      Inputs and Outputs
 * \defgroup macro   Macros
 * \defgroup math    Math
 * \defgroup md      MD Configuration
 * \defgroup nvidia  Nvidia GPU
 * \defgroup spatial Spatial Decomposition
 * \defgroup util    Utilities
 */


/**
 * \defgroup geom  Geometrical Restraints
 * \ingroup ff
 * \defgroup mpole Multipole Electrostatics
 * \ingroup ff
 * \defgroup vdw   Van der Waals (VDW)
 * \ingroup ff
 */


/**
 * \defgroup atomic Atomic Operations
 * \ingroup util
 * \defgroup ebuf   Energy Buffer
 * \ingroup util
 * \defgroup mem    Memory and Pointer
 * \ingroup util
 */


#include "macro.h"
TINKER_NAMESPACE_BEGIN
/**
 * \brief Wrappers for the compiler builtins.
 */
namespace builtin {}


/**
 * \brief Integer flags for different kinds of atom data and energy routines.
 */
namespace calc {}


/**
 * \brief Compile-time math functions. Some of them are also run-time functions.
 */
namespace ct {}


/**
 * \brief Math functions running in parallel.
 */
namespace parallel {}


/**
 * \brief Starting from 0 (which must be the `(0,0,0)` box), number the spatial
 * decomposition boxes as if they are stored in the `c[nx][ny][nz]` (or
 * `f(nz,ny,nx)`) array.
 */
namespace spatial_v1 {}


/**
 * \brief Starting from 0 (which must be the `(0,0,0)` box), number the spatial
 * decomposition boxes in a "binary tree" way.
 */
namespace spatial_v2 {}
TINKER_NAMESPACE_END
