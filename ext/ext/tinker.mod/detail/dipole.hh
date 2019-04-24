#ifndef TINKER_MOD_DIPOLE_HH_
#define TINKER_MOD_DIPOLE_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace dipole {
extern int& ndipole;
extern int*& idpl;
extern double*& bdpl;
extern double*& sdpl;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(dipole, ndipole);
extern "C" int* m_tinker_mod(dipole, idpl);
extern "C" double* m_tinker_mod(dipole, bdpl);
extern "C" double* m_tinker_mod(dipole, sdpl);

int& ndipole = m_tinker_mod(dipole, ndipole);
int*& idpl = m_tinker_mod(dipole, idpl);
double*& bdpl = m_tinker_mod(dipole, bdpl);
double*& sdpl = m_tinker_mod(dipole, sdpl);
#endif

} TINKER_NAMESPACE_END

#endif
