#ifndef TINKER_MOD_VIRIAL_HH_
#define TINKER_MOD_VIRIAL_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace virial {
extern double (&vir)[3][3];
extern int& use_virial;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(virial, vir)[3][3];
extern "C" int m_tinker_mod(virial, use_virial);

double (&vir)[3][3] = m_tinker_mod(virial, vir);
int& use_virial = m_tinker_mod(virial, use_virial);
#endif

} TINKER_NAMESPACE_END

#endif
