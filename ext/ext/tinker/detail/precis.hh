#ifndef TINKER_MOD_PRECIS_HH_
#define TINKER_MOD_PRECIS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace precis {
extern double& tiny;
extern double& small;
extern double& huge;

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(precis, tiny);
extern "C" double TINKER_MOD(precis, small);
extern "C" double TINKER_MOD(precis, huge);

double& tiny = TINKER_MOD(precis, tiny);
double& small = TINKER_MOD(precis, small);
double& huge = TINKER_MOD(precis, huge);
#endif
} TINKER_NAMESPACE_END

#endif
