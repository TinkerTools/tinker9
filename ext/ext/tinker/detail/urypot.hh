#pragma once

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace urypot {
extern double& cury;
extern double& qury;
extern double& ureyunit;

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(urypot, cury);
extern "C" double TINKER_MOD(urypot, qury);
extern "C" double TINKER_MOD(urypot, ureyunit);

double& cury = TINKER_MOD(urypot, cury);
double& qury = TINKER_MOD(urypot, qury);
double& ureyunit = TINKER_MOD(urypot, ureyunit);
#endif
} TINKER_NAMESPACE_END
