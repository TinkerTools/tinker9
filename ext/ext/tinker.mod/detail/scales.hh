#ifndef TINKER_MOD_SCALES_HH_
#define TINKER_MOD_SCALES_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace scales {
extern double*& scale;
extern int& set_scale;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(scales, scale);
extern "C" int m_tinker_mod(scales, set_scale);

double*& scale = m_tinker_mod(scales, scale);
int& set_scale = m_tinker_mod(scales, set_scale);
#endif

} TINKER_NAMESPACE_END

#endif
