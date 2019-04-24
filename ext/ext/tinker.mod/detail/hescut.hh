#ifndef TINKER_MOD_HESCUT_HH_
#define TINKER_MOD_HESCUT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace hescut {
extern double& hesscut;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(hescut, hesscut);

double& hesscut = m_tinker_mod(hescut, hesscut);
#endif

} TINKER_NAMESPACE_END

#endif
