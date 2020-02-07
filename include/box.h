#pragma once
#include "rc_man.h"
#include "realn.h"


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup box
 * \brief Periodic boundary conditions (PBC).
 *
 * Based on the lengths and angles, three axes can be expressed as three
 * vectors
 * \f$ \mathbf{A}^t := (ax, ay, az) = (a,  0,  0),  \f$
 * \f$ \mathbf{B}^t := (bx, by, bz) = (bx, by, 0),  \f$
 * \f$ \mathbf{C}^t := (cx, cy, cz) = (cx, cy, cz); \f$
 *
 * or a matrix V
 * \f$ = (\mathbf{A}, \mathbf{B}, \mathbf{C})
 * = \begin{pmatrix}
 * a & bx & cx \\
 * 0 & by & cy \\
 * 0 & 0  & cz
 * \end{pmatrix}. \f$
 *
 * `lvec`: Real space lattice vectors as matrix rows, is the transpose of
 * matrix V.
 *    - `lvec1 = lvec(:,1) = (a, bx, cx)`.
 *    - `lvec2 = lvec(:,2) = (0, by, cy)`.
 *    - `lvec3 = lvec(:,3) = (0, 0, cz)`.
 *    - `lvec` \f$ = \begin{pmatrix}
 *       a  & 0  & 0 \\
 *       bx & by & 0 \\
 *       cx & cy & cz
 *       \end{pmatrix}. \f$
 *    - Monoclinic representation:
 *       - \f$ \begin{pmatrix}
 *         a  & 0 & 0 \\
 *         0  & b & 0 \\
 *         cx & 0 & cz
 *         \end{pmatrix}. \f$
 *       - alpha and gamma are 90 degrees.
 *    - Orthogonal Fortran representation:
 *       - \f$ \begin{pmatrix}
 *         a & 0 & 0 \\
 *         0 & b & 0 \\
 *         0 & 0 & c
 *         \end{pmatrix}. \f$
 *       - alpha, beta, and gamma are 90 degrees.
 *
 * `recip`: Reciprocal lattice vectors as matrix columns, is the inverse of
 * `lvec`.
 *    - Triclinic representation:
 *       - \f$ \begin{pmatrix}
 *         rax & 0   & 0 \\
 *         ray & rby & 0 \\
 *         raz & rbz & rcz
 *         \end{pmatrix}. \f$
 *    - Fractional coordinates: \f$ (fx, fy, fz) = (xr, yr, zr) \cdot recip. \f$
 *    - \f$ (xr, yr, zr) = (fx, fy, fz) \cdot recip^{-1}
 *                       = (fx, fy, fz) \cdot lvec. \f$
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
TINKER_EXTERN real3 recipa;
/**
 * \ingroup box
 * \brief Reciprocal lattice vector b.
 */
TINKER_EXTERN real3 recipb;
/**
 * \ingroup box
 * \brief Reciprocal lattice vector c.
 */
TINKER_EXTERN real3 recipc;


TINKER_EXTERN real3 lvec1, lvec2, lvec3;


#define TINKER_IMAGE_PARAMS                                                    \
   real3 lvec1, real3 lvec2, real3 lvec3, real3 recipa, real3 recipb,          \
      real3 recipc
#define TINKER_IMAGE_ARGS lvec1, lvec2, lvec3, recipa, recipb, recipc


void box_data(rc_op);
void copyout_box_data(const Box*);
/**
 * \ingroup box
 * \brief Get the volume of the PBC box.
 * \note This function may calculate the volume on-the-fly if there is not an
 * internal variable to save the volume, so the calculation must be cheap.
 * \note The volume is undefined for the non-PBC situation.
 */
real volbox();
TINKER_NAMESPACE_END
