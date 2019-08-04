#ifndef TINKER_MOD_RXNPOT_HH_
#define TINKER_MOD_RXNPOT_HH_

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace rxnpot {
extern int& rfterms;
extern double& rfsize;
extern double& rfbulkd;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(rxnpot, rfterms);
extern "C" double TINKER_MOD(rxnpot, rfsize);
extern "C" double TINKER_MOD(rxnpot, rfbulkd);

int& rfterms = TINKER_MOD(rxnpot, rfterms);
double& rfsize = TINKER_MOD(rxnpot, rfsize);
double& rfbulkd = TINKER_MOD(rxnpot, rfbulkd);
#endif
} TINKER_NAMESPACE_END

#endif
