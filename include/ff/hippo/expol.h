#pragma once
#include "ff/energybuffer.h"
#include "tool/rcman.h"

namespace tinker {
void expolData(RcOp);

void alterpol();

enum class ExpolScr
{
   NONE,
   S2U,
   S2,
   G
};
}
