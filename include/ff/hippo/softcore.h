#pragma once
#include "ff/energybuffer.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup hippovdw
/// \brief HIPPO soft-core flags.
enum class Hvdw : int
{
   DECOUPLE = 0,   ///< VDW lambda type: decouple.
   ANNIHILATE = 1, ///< VDW lambda type: annihilate.
};
}

namespace tinker {
TINKER_EXTERN Hvdw hvcouple;

TINKER_EXTERN real hvlam;

/// \ingroup hippovdw
/// \brief State weighting values (lambda) of all atoms for Repulsion 
/// and Dispersion potentials.
TINKER_EXTERN int* hmut;
}
