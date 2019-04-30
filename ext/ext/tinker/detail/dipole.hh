#ifndef TINKER_MOD_DIPOLE_HH_
#define TINKER_MOD_DIPOLE_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace dipole {
extern int& ndipole;
extern int*& idpl;
extern double*& bdpl;
extern double*& sdpl;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(dipole, ndipole);
extern "C" int* TINKER_MOD(dipole, idpl);
extern "C" double* TINKER_MOD(dipole, bdpl);
extern "C" double* TINKER_MOD(dipole, sdpl);

int& ndipole = TINKER_MOD(dipole, ndipole);
int*& idpl = TINKER_MOD(dipole, idpl);
double*& bdpl = TINKER_MOD(dipole, bdpl);
double*& sdpl = TINKER_MOD(dipole, sdpl);
#endif

} TINKER_NAMESPACE_END

#endif
