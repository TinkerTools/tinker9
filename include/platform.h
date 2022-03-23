#pragma once
#include "mod/platform.h"
#include "tool/rcman.h"

namespace tinker {
constexpr int UNSET_PLTFM = 0x000;
constexpr int ACC_PLTFM = 0x001;
constexpr int CU_PLTFM = 0x002;

void platformData(RcOp);
}
