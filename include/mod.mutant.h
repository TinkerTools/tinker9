#pragma once
#include "macro.h"
#include "mod.vdwpot.h"

namespace tinker {
TINKER_EXTERN evdw_t vcouple;
/**
 * \ingroup vdw
 * \brief Exponential factor for soft core buffered 14-7 potential.
 */
TINKER_EXTERN real scexp;
/**
 * \ingroup vdw
 * \brief Scale factor \f$ \alpha \f$ for soft core buffered 14-7 potential.
 */
TINKER_EXTERN real scalpha;

//====================================================================//

TINKER_EXTERN real vlam;
/**
 * \ingroup vdw
 * \brief
 * State weighting values (lambda) of all atoms for van der Waals potentials.
 */
TINKER_EXTERN int* mut;
}
