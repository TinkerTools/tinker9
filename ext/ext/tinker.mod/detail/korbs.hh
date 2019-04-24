#ifndef TINKER_MOD_KORBS_HH_
#define TINKER_MOD_KORBS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace korbs {
const int maxnpi = 500;
const int maxnpi5 = 200;
const int maxnpi4 = 200;
extern double (&sslope)[maxnpi];
extern double (&sslope5)[maxnpi5];
extern double (&sslope4)[maxnpi4];
extern double (&tslope)[maxnpi];
extern double (&tslope5)[maxnpi5];
extern double (&tslope4)[maxnpi4];
extern double*& electron;
extern double*& ionize;
extern double*& repulse;
extern char (&kpi)[maxnpi][8];
extern char (&kpi5)[maxnpi5][8];
extern char (&kpi4)[maxnpi4][8];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(korbs, sslope)[maxnpi];
extern "C" double m_tinker_mod(korbs, sslope5)[maxnpi5];
extern "C" double m_tinker_mod(korbs, sslope4)[maxnpi4];
extern "C" double m_tinker_mod(korbs, tslope)[maxnpi];
extern "C" double m_tinker_mod(korbs, tslope5)[maxnpi5];
extern "C" double m_tinker_mod(korbs, tslope4)[maxnpi4];
extern "C" double* m_tinker_mod(korbs, electron);
extern "C" double* m_tinker_mod(korbs, ionize);
extern "C" double* m_tinker_mod(korbs, repulse);
extern "C" char m_tinker_mod(korbs, kpi)[maxnpi][8];
extern "C" char m_tinker_mod(korbs, kpi5)[maxnpi5][8];
extern "C" char m_tinker_mod(korbs, kpi4)[maxnpi4][8];

double (&sslope)[maxnpi] = m_tinker_mod(korbs, sslope);
double (&sslope5)[maxnpi5] = m_tinker_mod(korbs, sslope5);
double (&sslope4)[maxnpi4] = m_tinker_mod(korbs, sslope4);
double (&tslope)[maxnpi] = m_tinker_mod(korbs, tslope);
double (&tslope5)[maxnpi5] = m_tinker_mod(korbs, tslope5);
double (&tslope4)[maxnpi4] = m_tinker_mod(korbs, tslope4);
double*& electron = m_tinker_mod(korbs, electron);
double*& ionize = m_tinker_mod(korbs, ionize);
double*& repulse = m_tinker_mod(korbs, repulse);
char (&kpi)[maxnpi][8] = m_tinker_mod(korbs, kpi);
char (&kpi5)[maxnpi5][8] = m_tinker_mod(korbs, kpi5);
char (&kpi4)[maxnpi4][8] = m_tinker_mod(korbs, kpi4);
#endif

} TINKER_NAMESPACE_END

#endif
