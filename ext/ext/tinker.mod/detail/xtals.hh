#ifndef TINKER_MOD_XTALS_HH_
#define TINKER_MOD_XTALS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace xtals {
const int maxlsq = 1000;
const int maxrsd = 1000;
extern int& nxtal;
extern int& nvary;
extern int (&ivary)[maxlsq];
extern int (&iresid)[maxrsd];
extern int (&vary)[maxlsq][2];
extern double& e0_lattice;
extern char (&vartyp)[maxlsq][16];
extern char (&rsdtyp)[maxrsd][16];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(xtals, nxtal);
extern "C" int m_tinker_mod(xtals, nvary);
extern "C" int m_tinker_mod(xtals, ivary)[maxlsq];
extern "C" int m_tinker_mod(xtals, iresid)[maxrsd];
extern "C" int m_tinker_mod(xtals, vary)[maxlsq][2];
extern "C" double m_tinker_mod(xtals, e0_lattice);
extern "C" char m_tinker_mod(xtals, vartyp)[maxlsq][16];
extern "C" char m_tinker_mod(xtals, rsdtyp)[maxrsd][16];

int& nxtal = m_tinker_mod(xtals, nxtal);
int& nvary = m_tinker_mod(xtals, nvary);
int (&ivary)[maxlsq] = m_tinker_mod(xtals, ivary);
int (&iresid)[maxrsd] = m_tinker_mod(xtals, iresid);
int (&vary)[maxlsq][2] = m_tinker_mod(xtals, vary);
double& e0_lattice = m_tinker_mod(xtals, e0_lattice);
char (&vartyp)[maxlsq][16] = m_tinker_mod(xtals, vartyp);
char (&rsdtyp)[maxrsd][16] = m_tinker_mod(xtals, rsdtyp);
#endif

} TINKER_NAMESPACE_END

#endif
