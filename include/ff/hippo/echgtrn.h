#pragma once
#include "tool/rcman.h"

namespace tinker {
/// \ingroup chgtrn
enum class Chgtrn
{
   SEPARATE,
   COMBINED
};

/// \ingroup chgtrn
void echgtrnData(RcOp);
/// \ingroup chgtrn
void echgtrn(int vers);
}
