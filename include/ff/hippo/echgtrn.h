#pragma once
#include "tool/rcman.h"

namespace tinker {
enum class Chgtrn
{
   SEPARATE,
   COMBINED
};

void echgtrnData(RcOp);
void echgtrn(int vers);
}
