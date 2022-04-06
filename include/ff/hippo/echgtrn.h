#pragma once
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
