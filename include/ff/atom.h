#pragma once
#include "precision.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup ff
void nData(RcOp);
/// \ingroup ff
void massData(RcOp);
/// \ingroup ff
void xyzData(RcOp);

/// \ingroup ff
/// \brief Update #x, #y, #z by #xpos, #ypos, and #zpos.
/// If #xpos etc. are only aliases, return directly.
void copyPosToXyz();

/// \ingroup ff
/// \brief Update #x, #y, #z by #xpos, #ypos, and #zpos.
/// If #xpos etc. are only aliases, return directly.
/// \param refreshNBList  If `true`, refresh the neighbor lists by the end.
void copyPosToXyz(bool refreshNBList);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
/// \ingroup ff
/// \brief Number of atoms.
TINKER_EXTERN int padded_n;
/// \ingroup ff
/// \brief Number of atoms padded by #WARP_SIZE.
/// \see WARP_SIZE
TINKER_EXTERN int n;
/// \ingroup ff
/// \brief Number of the trajectory frames.
TINKER_EXTERN int trajn;

/// \ingroup ff
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

/// \ingroup ff
/// \brief Atomic mass.
TINKER_EXTERN double* mass;
/// \ingroup ff
/// \brief Inversed atomic mass.
TINKER_EXTERN double* massinv;
}
