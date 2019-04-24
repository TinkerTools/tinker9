#ifndef TINKER_MOD_RXNPOT_HH_
#define TINKER_MOD_RXNPOT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace rxnpot {
extern int& rfterms;
extern double& rfsize;
extern double& rfbulkd;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(rxnpot, rfterms);
extern "C" double m_tinker_mod(rxnpot, rfsize);
extern "C" double m_tinker_mod(rxnpot, rfbulkd);

int& rfterms = m_tinker_mod(rxnpot, rfterms);
double& rfsize = m_tinker_mod(rxnpot, rfsize);
double& rfbulkd = m_tinker_mod(rxnpot, rfbulkd);
#endif

} TINKER_NAMESPACE_END

#endif
