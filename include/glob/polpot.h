#pragma once
#include "macro.h"

// polpot
namespace tinker {
TINKER_EXTERN real u1scale;
TINKER_EXTERN real u2scale;
TINKER_EXTERN real u3scale;
TINKER_EXTERN real u4scale;
TINKER_EXTERN real d1scale;
TINKER_EXTERN real d2scale;
TINKER_EXTERN real d3scale;
TINKER_EXTERN real d4scale;
TINKER_EXTERN real p2scale;
TINKER_EXTERN real p3scale;
TINKER_EXTERN real p4scale;
TINKER_EXTERN real p5scale;
TINKER_EXTERN real p2iscale;
TINKER_EXTERN real p3iscale;
TINKER_EXTERN real p4iscale;
TINKER_EXTERN real p5iscale;

TINKER_EXTERN real udiag;

//====================================================================//

TINKER_EXTERN int nuexclude;
TINKER_EXTERN int (*uexclude)[2];
TINKER_EXTERN real* uexclude_scale;

TINKER_EXTERN int ndpexclude;
TINKER_EXTERN int (*dpexclude)[2];
TINKER_EXTERN real (*dpexclude_scale)[2];

TINKER_EXTERN int ndpuexclude;
TINKER_EXTERN int (*dpuexclude)[2];
TINKER_EXTERN real (*dpuexclude_scale)[3];
}
