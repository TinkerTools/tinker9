#ifndef TINKER_MOD_BNDPOT_HH_
#define TINKER_MOD_BNDPOT_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace bndpot {
extern double& cbnd;
extern double& qbnd;
extern double& bndunit;
extern char (&bndtyp)[8];

#ifdef TINKER_MOD_CPP_
extern "C" double TINKER_MOD(bndpot, cbnd);
extern "C" double TINKER_MOD(bndpot, qbnd);
extern "C" double TINKER_MOD(bndpot, bndunit);
extern "C" char TINKER_MOD(bndpot, bndtyp)[8];

double& cbnd = TINKER_MOD(bndpot, cbnd);
double& qbnd = TINKER_MOD(bndpot, qbnd);
double& bndunit = TINKER_MOD(bndpot, bndunit);
char (&bndtyp)[8] = TINKER_MOD(bndpot, bndtyp);
#endif
} TINKER_NAMESPACE_END

#endif
