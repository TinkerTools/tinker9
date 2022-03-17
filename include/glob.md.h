#pragma once
#include "precision.h"

namespace tinker {
/// \ingroup mdcalc
/// \var rc_flag
/// \brief Global bitmask.
TINKER_EXTERN int rc_flag;

/// \ingroup mdpq
/// \{
/// \var n
/// \brief Number of atoms.
///
/// \var padded_n
/// \brief Number of atoms padded by #WARP_SIZE.
/// \see WARP_SIZE
///
/// \var trajn
/// \brief Number of the trajectory frames.
/// \}
TINKER_EXTERN int n;
TINKER_EXTERN int padded_n;
TINKER_EXTERN int trajn;

/// \ingroup mdpq
/// \{
/// \var x
/// \brief Current coordinates used in energy evaluation and neighbor lists.
/// \var y
/// \copydoc x
/// \var z
/// \copydoc x
///
/// \var trajx
/// \brief Coordinates of all the trajectory frames.
/// \var trajy
/// \copydoc trajx
/// \var trajz
/// \copydoc trajx
///
/// \var xpos
/// \brief Coordinates used in integrators.
/// \note
///    - New arrays will be allocated only if `sizeof(pos_prec) > sizeof(real)`,
///    otherwise, they will be aliases of #x, #y, and #z.
///    - Whenever #xpos, #ypos, #zpos get updated by integrators, barostats etc.,
///    #x, #y, #z must be updated immediately.
/// \see pos_prec
/// \see real
/// \var ypos
/// \copydoc xpos
/// \var zpos
/// \copydoc xpos
/// \}
TINKER_EXTERN real *x, *y, *z;
TINKER_EXTERN real *trajx, *trajy, *trajz;
TINKER_EXTERN pos_prec *xpos, *ypos, *zpos;
static_assert(sizeof(pos_prec) >= sizeof(real), "Type pos_prec cannot be shorter than type real.");

/// \ingroup mdpq
/// \{
/// \var mass
/// \brief Atomic mass.
/// \var massinv
/// \brief Inversed atomic mass.
///
/// \var vx
/// \brief Velocities.
/// \var vy
/// \copydoc vx
/// \var vz
/// \copydoc vx
/// \}
TINKER_EXTERN double* mass;
TINKER_EXTERN double* massinv;
TINKER_EXTERN vel_prec *vx, *vy, *vz;

/// \ingroup mdintg
/// \{
/// \var gx1
/// \brief Gradient for the fast RESPA energy terms.
/// \var gy1
/// \copydoc gx1
/// \var gz1
/// \copydoc gx1
/// \var gx2
/// \brief Gradient for the slow RESPA energy terms.
/// \var gy2
/// \copydoc gx2
/// \var gz2
/// \copydoc gx2
/// \}
TINKER_EXTERN grad_prec *gx1, *gy1, *gz1;
TINKER_EXTERN grad_prec *gx2, *gy2, *gz2;

/// \ingroup mdintg
/// \{
/// \var x_pmonte
/// \brief Temporary coordinates created for the Monte Carlo barostat.
/// \var y_pmonte
/// \copydoc x_pmonte
/// \var z_pmonte
/// \copydoc x_pmonte
/// \}
TINKER_EXTERN pos_prec *x_pmonte, *y_pmonte, *z_pmonte;
}
