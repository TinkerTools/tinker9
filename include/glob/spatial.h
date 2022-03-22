#pragma once
#include "ff/spatial.h"
#include "ff/spatial2.h"
#include "macro.h"

namespace tinker {
TINKER_EXTERN Spatial2Unit cspatial_v2_unit;
TINKER_EXTERN Spatial2Unit vspatial_v2_unit;
TINKER_EXTERN Spatial2Unit uspatial_v2_unit;
TINKER_EXTERN Spatial2Unit mspatial_v2_unit;
TINKER_EXTERN Spatial2Unit dspspatial_v2_unit;

constexpr int cspatial_fresh_mask_echglj = 0x00000001;
}
