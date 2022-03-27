#pragma once
#include "ff/nblist.h"
#include "ff/spatial.h"

// nblist
namespace tinker {
TINKER_EXTERN NBListUnit vlist_unit;
TINKER_EXTERN NBListUnit clist_unit;
TINKER_EXTERN NBListUnit mlist_unit;
TINKER_EXTERN NBListUnit ulist_unit;
TINKER_EXTERN NBListUnit dsplist_unit;
}

// spatial
namespace tinker {
TINKER_EXTERN SpatialUnit cspatial_v2_unit;
TINKER_EXTERN SpatialUnit vspatial_v2_unit;
TINKER_EXTERN SpatialUnit uspatial_v2_unit;
TINKER_EXTERN SpatialUnit mspatial_v2_unit;
TINKER_EXTERN SpatialUnit dspspatial_v2_unit;

constexpr int cspatial_fresh_mask_echglj = 0x00000001;
}
