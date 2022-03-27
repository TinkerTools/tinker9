#pragma once
#include "precision.h"

namespace tinker {
/// \ingroup rc
/// \var rc_flag
/// \brief Global bitmask.
TINKER_EXTERN int rc_flag;

/// \ingroup mdpq
/// \{
/// \var vx
/// \brief Velocities.
/// \var vy
/// \copydoc vx
/// \var vz
/// \copydoc vx
/// \}
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

TINKER_EXTERN double qbar;
TINKER_EXTERN double vbar;
TINKER_EXTERN double vbar_matrix[3][3];
}
