#ifndef TINKER_MOD_ATOMS_HH_
#define TINKER_MOD_ATOMS_HH_

#include "util/macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace atoms {
using namespace sizes;

extern int& n;
extern int (&type)[maxatm];
extern double (&x)[maxatm];
extern double (&y)[maxatm];
extern double (&z)[maxatm];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(atoms, n);
extern "C" int m_tinker_mod(atoms, type)[maxatm];
extern "C" double m_tinker_mod(atoms, x)[maxatm];
extern "C" double m_tinker_mod(atoms, y)[maxatm];
extern "C" double m_tinker_mod(atoms, z)[maxatm];

int& n = m_tinker_mod(atoms, n);
int (&type)[maxatm] = m_tinker_mod(atoms, type);
double (&x)[maxatm] = m_tinker_mod(atoms, x);
double (&y)[maxatm] = m_tinker_mod(atoms, y);
double (&z)[maxatm] = m_tinker_mod(atoms, z);
#endif

} TINKER_NAMESPACE_END

#endif
