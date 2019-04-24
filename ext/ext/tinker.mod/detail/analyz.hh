#ifndef TINKER_MOD_ANALYZ_HH_
#define TINKER_MOD_ANALYZ_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace analyz {
extern double*& aesum;
extern double*& aeb;
extern double*& aea;
extern double*& aeba;
extern double*& aeub;
extern double*& aeaa;
extern double*& aeopb;
extern double*& aeopd;
extern double*& aeid;
extern double*& aeit;
extern double*& aet;
extern double*& aept;
extern double*& aebt;
extern double*& aeat;
extern double*& aett;
extern double*& aev;
extern double*& aer;
extern double*& aedsp;
extern double*& aec;
extern double*& aecd;
extern double*& aed;
extern double*& aem;
extern double*& aep;
extern double*& aect;
extern double*& aerxf;
extern double*& aes;
extern double*& aelf;
extern double*& aeg;
extern double*& aex;

#ifdef TINKER_MOD_CPP_
extern "C" double* m_tinker_mod(analyz, aesum);
extern "C" double* m_tinker_mod(analyz, aeb);
extern "C" double* m_tinker_mod(analyz, aea);
extern "C" double* m_tinker_mod(analyz, aeba);
extern "C" double* m_tinker_mod(analyz, aeub);
extern "C" double* m_tinker_mod(analyz, aeaa);
extern "C" double* m_tinker_mod(analyz, aeopb);
extern "C" double* m_tinker_mod(analyz, aeopd);
extern "C" double* m_tinker_mod(analyz, aeid);
extern "C" double* m_tinker_mod(analyz, aeit);
extern "C" double* m_tinker_mod(analyz, aet);
extern "C" double* m_tinker_mod(analyz, aept);
extern "C" double* m_tinker_mod(analyz, aebt);
extern "C" double* m_tinker_mod(analyz, aeat);
extern "C" double* m_tinker_mod(analyz, aett);
extern "C" double* m_tinker_mod(analyz, aev);
extern "C" double* m_tinker_mod(analyz, aer);
extern "C" double* m_tinker_mod(analyz, aedsp);
extern "C" double* m_tinker_mod(analyz, aec);
extern "C" double* m_tinker_mod(analyz, aecd);
extern "C" double* m_tinker_mod(analyz, aed);
extern "C" double* m_tinker_mod(analyz, aem);
extern "C" double* m_tinker_mod(analyz, aep);
extern "C" double* m_tinker_mod(analyz, aect);
extern "C" double* m_tinker_mod(analyz, aerxf);
extern "C" double* m_tinker_mod(analyz, aes);
extern "C" double* m_tinker_mod(analyz, aelf);
extern "C" double* m_tinker_mod(analyz, aeg);
extern "C" double* m_tinker_mod(analyz, aex);

double*& aesum = m_tinker_mod(analyz, aesum);
double*& aeb = m_tinker_mod(analyz, aeb);
double*& aea = m_tinker_mod(analyz, aea);
double*& aeba = m_tinker_mod(analyz, aeba);
double*& aeub = m_tinker_mod(analyz, aeub);
double*& aeaa = m_tinker_mod(analyz, aeaa);
double*& aeopb = m_tinker_mod(analyz, aeopb);
double*& aeopd = m_tinker_mod(analyz, aeopd);
double*& aeid = m_tinker_mod(analyz, aeid);
double*& aeit = m_tinker_mod(analyz, aeit);
double*& aet = m_tinker_mod(analyz, aet);
double*& aept = m_tinker_mod(analyz, aept);
double*& aebt = m_tinker_mod(analyz, aebt);
double*& aeat = m_tinker_mod(analyz, aeat);
double*& aett = m_tinker_mod(analyz, aett);
double*& aev = m_tinker_mod(analyz, aev);
double*& aer = m_tinker_mod(analyz, aer);
double*& aedsp = m_tinker_mod(analyz, aedsp);
double*& aec = m_tinker_mod(analyz, aec);
double*& aecd = m_tinker_mod(analyz, aecd);
double*& aed = m_tinker_mod(analyz, aed);
double*& aem = m_tinker_mod(analyz, aem);
double*& aep = m_tinker_mod(analyz, aep);
double*& aect = m_tinker_mod(analyz, aect);
double*& aerxf = m_tinker_mod(analyz, aerxf);
double*& aes = m_tinker_mod(analyz, aes);
double*& aelf = m_tinker_mod(analyz, aelf);
double*& aeg = m_tinker_mod(analyz, aeg);
double*& aex = m_tinker_mod(analyz, aex);
#endif

} TINKER_NAMESPACE_END

#endif
