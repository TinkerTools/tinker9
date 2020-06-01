#pragma once
#include "macro.h"


namespace tinker {
/// \ingroup geom
/// \brief Number of group distance restraints to be applied.
TINKER_EXTERN int ngfix;
/// \ingroup geom
/// \brief Group numbers defining each group distance restraint.
TINKER_EXTERN int (*igfix)[2];
/// \ingroup geom
/// \brief Force constant and target range for each group distance.
TINKER_EXTERN real (*gfix)[3];
}
