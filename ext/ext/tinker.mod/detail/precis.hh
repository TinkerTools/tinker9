#ifndef TINKER_MOD_PRECIS_HH_
#define TINKER_MOD_PRECIS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace precis {
extern double& tiny;
extern double& small;
extern double& huge;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(precis, tiny);
extern "C" double m_tinker_mod(precis, small);
extern "C" double m_tinker_mod(precis, huge);

double& tiny = m_tinker_mod(precis, tiny);
double& small = m_tinker_mod(precis, small);
double& huge = m_tinker_mod(precis, huge);
#endif

} TINKER_NAMESPACE_END

#endif
