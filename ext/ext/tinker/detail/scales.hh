#pragma once

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace scales {
extern double*& scale;
extern int& set_scale;

#ifdef TINKER_MOD_CPP_
extern "C" double* TINKER_MOD(scales, scale);
extern "C" int TINKER_MOD(scales, set_scale);

double*& scale = TINKER_MOD(scales, scale);
int& set_scale = TINKER_MOD(scales, set_scale);
#endif
} TINKER_NAMESPACE_END
