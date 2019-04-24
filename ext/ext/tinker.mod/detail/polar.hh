#ifndef TINKER_MOD_POLAR_HH_
#define TINKER_MOD_POLAR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace polar {
extern int& npolar;
extern int*& ipolar;
extern double*& polarity;
extern double*& thole;
extern double*& pdamp;
extern double*& udir;
extern double*& udirp;
extern double*& udirs;
extern double*& udirps;
extern double*& uind;
extern double*& uinp;
extern double*& uinds;
extern double*& uinps;
extern double*& uexact;
extern int*& douind;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(polar, npolar);
extern "C" int* m_tinker_mod(polar, ipolar);
extern "C" double* m_tinker_mod(polar, polarity);
extern "C" double* m_tinker_mod(polar, thole);
extern "C" double* m_tinker_mod(polar, pdamp);
extern "C" double* m_tinker_mod(polar, udir);
extern "C" double* m_tinker_mod(polar, udirp);
extern "C" double* m_tinker_mod(polar, udirs);
extern "C" double* m_tinker_mod(polar, udirps);
extern "C" double* m_tinker_mod(polar, uind);
extern "C" double* m_tinker_mod(polar, uinp);
extern "C" double* m_tinker_mod(polar, uinds);
extern "C" double* m_tinker_mod(polar, uinps);
extern "C" double* m_tinker_mod(polar, uexact);
extern "C" int* m_tinker_mod(polar, douind);

int& npolar = m_tinker_mod(polar, npolar);
int*& ipolar = m_tinker_mod(polar, ipolar);
double*& polarity = m_tinker_mod(polar, polarity);
double*& thole = m_tinker_mod(polar, thole);
double*& pdamp = m_tinker_mod(polar, pdamp);
double*& udir = m_tinker_mod(polar, udir);
double*& udirp = m_tinker_mod(polar, udirp);
double*& udirs = m_tinker_mod(polar, udirs);
double*& udirps = m_tinker_mod(polar, udirps);
double*& uind = m_tinker_mod(polar, uind);
double*& uinp = m_tinker_mod(polar, uinp);
double*& uinds = m_tinker_mod(polar, uinds);
double*& uinps = m_tinker_mod(polar, uinps);
double*& uexact = m_tinker_mod(polar, uexact);
int*& douind = m_tinker_mod(polar, douind);
#endif

} TINKER_NAMESPACE_END

#endif
