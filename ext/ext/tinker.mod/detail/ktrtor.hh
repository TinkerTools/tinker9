#ifndef TINKER_MOD_KTRTOR_HH_
#define TINKER_MOD_KTRTOR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace ktrtor {
const int maxntt = 100;
const int maxtgrd = 30;
const int maxtgrd2 = maxtgrd*maxtgrd;
extern int (&tnx)[maxntt];
extern int (&tny)[maxntt];
extern double (&ttx)[maxntt][maxtgrd];
extern double (&tty)[maxntt][maxtgrd];
extern double (&tbf)[maxntt][maxtgrd2];
extern double (&tbx)[maxntt][maxtgrd2];
extern double (&tby)[maxntt][maxtgrd2];
extern double (&tbxy)[maxntt][maxtgrd2];
extern char (&ktt)[maxntt][20];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(ktrtor, tnx)[maxntt];
extern "C" int m_tinker_mod(ktrtor, tny)[maxntt];
extern "C" double m_tinker_mod(ktrtor, ttx)[maxntt][maxtgrd];
extern "C" double m_tinker_mod(ktrtor, tty)[maxntt][maxtgrd];
extern "C" double m_tinker_mod(ktrtor, tbf)[maxntt][maxtgrd2];
extern "C" double m_tinker_mod(ktrtor, tbx)[maxntt][maxtgrd2];
extern "C" double m_tinker_mod(ktrtor, tby)[maxntt][maxtgrd2];
extern "C" double m_tinker_mod(ktrtor, tbxy)[maxntt][maxtgrd2];
extern "C" char m_tinker_mod(ktrtor, ktt)[maxntt][20];

int (&tnx)[maxntt] = m_tinker_mod(ktrtor, tnx);
int (&tny)[maxntt] = m_tinker_mod(ktrtor, tny);
double (&ttx)[maxntt][maxtgrd] = m_tinker_mod(ktrtor, ttx);
double (&tty)[maxntt][maxtgrd] = m_tinker_mod(ktrtor, tty);
double (&tbf)[maxntt][maxtgrd2] = m_tinker_mod(ktrtor, tbf);
double (&tbx)[maxntt][maxtgrd2] = m_tinker_mod(ktrtor, tbx);
double (&tby)[maxntt][maxtgrd2] = m_tinker_mod(ktrtor, tby);
double (&tbxy)[maxntt][maxtgrd2] = m_tinker_mod(ktrtor, tbxy);
char (&ktt)[maxntt][20] = m_tinker_mod(ktrtor, ktt);
#endif

} TINKER_NAMESPACE_END

#endif
