#ifndef TINKER_MOD_MPOLE_HH_
#define TINKER_MOD_MPOLE_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace mpole {
const int maxpole = 13;
extern int& npole;
extern int*& ipole;
extern int*& polsiz;
extern int*& pollist;
extern int*& zaxis;
extern int*& xaxis;
extern int*& yaxis;
extern double*& pole;
extern double*& rpole;
extern double*& spole;
extern double*& srpole;
extern char (*&polaxe)[8];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(mpole, npole);
extern "C" int* m_tinker_mod(mpole, ipole);
extern "C" int* m_tinker_mod(mpole, polsiz);
extern "C" int* m_tinker_mod(mpole, pollist);
extern "C" int* m_tinker_mod(mpole, zaxis);
extern "C" int* m_tinker_mod(mpole, xaxis);
extern "C" int* m_tinker_mod(mpole, yaxis);
extern "C" double* m_tinker_mod(mpole, pole);
extern "C" double* m_tinker_mod(mpole, rpole);
extern "C" double* m_tinker_mod(mpole, spole);
extern "C" double* m_tinker_mod(mpole, srpole);
extern "C" char (*m_tinker_mod(mpole, polaxe))[8];

int& npole = m_tinker_mod(mpole, npole);
int*& ipole = m_tinker_mod(mpole, ipole);
int*& polsiz = m_tinker_mod(mpole, polsiz);
int*& pollist = m_tinker_mod(mpole, pollist);
int*& zaxis = m_tinker_mod(mpole, zaxis);
int*& xaxis = m_tinker_mod(mpole, xaxis);
int*& yaxis = m_tinker_mod(mpole, yaxis);
double*& pole = m_tinker_mod(mpole, pole);
double*& rpole = m_tinker_mod(mpole, rpole);
double*& spole = m_tinker_mod(mpole, spole);
double*& srpole = m_tinker_mod(mpole, srpole);
char (*&polaxe)[8] = m_tinker_mod(mpole, polaxe);
#endif

} TINKER_NAMESPACE_END

#endif
