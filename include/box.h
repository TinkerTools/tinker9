#pragma once
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup box
 * \brief Periodic boundary conditions (PBC).
 *
 * `lvec[3][3]`: Real space lattice vectors as matrix rows.
 *    - Fortran representation:
 *       - \f$\begin{pmatrix}
 *         a1 & a2 & a3 \\
 *         b1 & b2 & b3 \\
 *         c1 & c2 & c3
 *         \end{pmatrix}\f$
 *    - C representation:
 *       - \code
 *         {{a1, b1, c1},
 *          {a2, b2, c2},
 *          {a3, b3, c3}};
 *         \endcode
 *    - Triclinic Fortran representation:
 *       - \f$\begin{pmatrix}
 *         a1 & 0  & 0 \\
 *         b1 & b2 & 0 \\
 *         c1 & c2 & c3
 *         \end{pmatrix}\f$
 *    - Monoclinic Fortran representation:
 *       - \f$\begin{pmatrix}
 *         a1 & 0  & 0 \\
 *         0  & b2 & 0 \\
 *         c1 & 0  & c3
 *         \end{pmatrix}\f$
 *       - alpha and gamma are 90 degrees.
 *    - Orthogonal Fortran representation:
 *       - \f$\begin{pmatrix}
 *         a1 & 0  & 0 \\
 *         0  & b2 & 0 \\
 *         0  & 0  & c3
 *         \end{pmatrix}\f$
 *       - alpha, beta, and gamma are 90 degrees.
 *
 * `recip[3][3]`: Reciprocal lattice vectors as matrix columns;
 * the inverse of `lev[3][3]`.
 *    - Fortran representation:
 *       - \f$\begin{pmatrix}
 *         ra1 & rb1 & rc1 \\
 *         ra2 & rb2 & rc2 \\
 *         ra3 & rb3 & rc3
 *         \end{pmatrix}\f$
 *    - Fractional coordinates: \f$ f_{x,y,z} = r_{a,b,c} \cdot r \f$.
 */
struct Box
{
   /// \brief Shape of the PBC box.
   typedef enum
   {
      null = 0x000,  ///< Do not use PBC.
      ortho = 0x001, ///< Orthogonal PBC.
      mono = 0x002,  ///< Monoclinic PBC.
      tri = 0x004,   ///< Triclinic PBC.
      oct = 0x008    ///< Truncated octahedron PBC.
   } Shape;
   real lvec[3][3];
   real recip[3][3];
   real volbox; ///< Volume of the PBC box.
   Shape shape; ///< Shape of the PBC box
};

/**
 * \ingroup box
 * \brief Device pointer to the PBC box.
 */
TINKER_EXTERN Box* box;
/**
 * \ingroup box
 * \brief Device pointer to the current PBC box of a trajectory.
 */
TINKER_EXTERN Box* trajbox;
/**
 * \ingroup box
 * \brief Reciprocal lattice vector a.
 */
TINKER_EXTERN real3 recip_a;
/**
 * \ingroup box
 * \brief Reciprocal lattice vector b.
 */
TINKER_EXTERN real3 recip_b;
/**
 * \ingroup box
 * \brief Reciprocal lattice vector c.
 */
TINKER_EXTERN real3 recip_c;

void box_data(rc_op);
void copyout_box_data(const Box*);
TINKER_NAMESPACE_END
