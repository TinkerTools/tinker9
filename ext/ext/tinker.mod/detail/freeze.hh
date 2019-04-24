#ifndef TINKER_MOD_FREEZE_HH_
#define TINKER_MOD_FREEZE_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace freeze {
extern int& nrat;
extern int& nratx;
extern int*& iratx;
extern int*& kratx;
extern int*& irat;
extern double& rateps;
extern double*& krat;
extern int& use_rattle;
extern int*& ratimage;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(freeze, nrat);
extern "C" int m_tinker_mod(freeze, nratx);
extern "C" int* m_tinker_mod(freeze, iratx);
extern "C" int* m_tinker_mod(freeze, kratx);
extern "C" int* m_tinker_mod(freeze, irat);
extern "C" double m_tinker_mod(freeze, rateps);
extern "C" double* m_tinker_mod(freeze, krat);
extern "C" int m_tinker_mod(freeze, use_rattle);
extern "C" int* m_tinker_mod(freeze, ratimage);

int& nrat = m_tinker_mod(freeze, nrat);
int& nratx = m_tinker_mod(freeze, nratx);
int*& iratx = m_tinker_mod(freeze, iratx);
int*& kratx = m_tinker_mod(freeze, kratx);
int*& irat = m_tinker_mod(freeze, irat);
double& rateps = m_tinker_mod(freeze, rateps);
double*& krat = m_tinker_mod(freeze, krat);
int& use_rattle = m_tinker_mod(freeze, use_rattle);
int*& ratimage = m_tinker_mod(freeze, ratimage);
#endif

} TINKER_NAMESPACE_END

#endif
