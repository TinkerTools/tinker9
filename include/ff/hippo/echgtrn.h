#pragma once
#include "ff/evdw.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup chgtrn
enum class Chgtrn
{
   NONE,
   SEPARATE, // hippo
   COMBINED  // amoeba+
};

/// \ingroup chgtrn
void echgtrnData(RcOp);
/// \ingroup chgtrn
void echgtrn(int vers);
}

namespace tinker {
/// \ingroup chgtrn
/// \brief state weighting value for electrostatic potentials
TINKER_EXTERN real elam;
}
