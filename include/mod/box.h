#pragma once
#include "ff/box1.h"

namespace tinker {
TINKER_EXTERN BoxShape box_shape;
TINKER_EXTERN real3 lvec1;
TINKER_EXTERN real3 lvec2;
TINKER_EXTERN real3 lvec3;
TINKER_EXTERN real3 recipa;
TINKER_EXTERN real3 recipb;
TINKER_EXTERN real3 recipc;

#define TINKER_IMAGE_LVEC_PARAMS  real3 lvec1, real3 lvec2, real3 lvec3
#define TINKER_IMAGE_LVEC_ARGS    lvec1, lvec2, lvec3
#define TINKER_IMAGE_RECIP_PARAMS real3 recipa, real3 recipb, real3 recipc
#define TINKER_IMAGE_RECIP_ARGS   recipa, recipb, recipc
#define TINKER_IMAGE_PARAMS       BoxShape box_shape, TINKER_IMAGE_LVEC_PARAMS, TINKER_IMAGE_RECIP_PARAMS
#define TINKER_IMAGE_ARGS         box_shape, TINKER_IMAGE_LVEC_ARGS, TINKER_IMAGE_RECIP_ARGS

/// \ingroup box
/// \brief Host pointer to the PBC boxes of a trajectory.
TINKER_EXTERN Box* trajbox;
}
