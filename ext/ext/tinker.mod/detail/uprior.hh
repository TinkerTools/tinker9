#ifndef TINKER_MOD_UPRIOR_HH_
#define TINKER_MOD_UPRIOR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace uprior {
const int maxpred = 17;
extern int& nualt;
extern int& maxualt;
extern double (&gear)[maxpred];
extern double (&aspc)[maxpred];
extern double (&bpred)[maxpred];
extern double (&bpredp)[maxpred];
extern double (&bpreds)[maxpred];
extern double (&bpredps)[maxpred];
extern double*& udalt;
extern double*& upalt;
extern double*& usalt;
extern double*& upsalt;
extern int& use_pred;
extern char (&polpred)[4];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(uprior, nualt);
extern "C" int m_tinker_mod(uprior, maxualt);
extern "C" double m_tinker_mod(uprior, gear)[maxpred];
extern "C" double m_tinker_mod(uprior, aspc)[maxpred];
extern "C" double m_tinker_mod(uprior, bpred)[maxpred];
extern "C" double m_tinker_mod(uprior, bpredp)[maxpred];
extern "C" double m_tinker_mod(uprior, bpreds)[maxpred];
extern "C" double m_tinker_mod(uprior, bpredps)[maxpred];
extern "C" double* m_tinker_mod(uprior, udalt);
extern "C" double* m_tinker_mod(uprior, upalt);
extern "C" double* m_tinker_mod(uprior, usalt);
extern "C" double* m_tinker_mod(uprior, upsalt);
extern "C" int m_tinker_mod(uprior, use_pred);
extern "C" char m_tinker_mod(uprior, polpred)[4];

int& nualt = m_tinker_mod(uprior, nualt);
int& maxualt = m_tinker_mod(uprior, maxualt);
double (&gear)[maxpred] = m_tinker_mod(uprior, gear);
double (&aspc)[maxpred] = m_tinker_mod(uprior, aspc);
double (&bpred)[maxpred] = m_tinker_mod(uprior, bpred);
double (&bpredp)[maxpred] = m_tinker_mod(uprior, bpredp);
double (&bpreds)[maxpred] = m_tinker_mod(uprior, bpreds);
double (&bpredps)[maxpred] = m_tinker_mod(uprior, bpredps);
double*& udalt = m_tinker_mod(uprior, udalt);
double*& upalt = m_tinker_mod(uprior, upalt);
double*& usalt = m_tinker_mod(uprior, usalt);
double*& upsalt = m_tinker_mod(uprior, upsalt);
int& use_pred = m_tinker_mod(uprior, use_pred);
char (&polpred)[4] = m_tinker_mod(uprior, polpred);
#endif

} TINKER_NAMESPACE_END

#endif
