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
extern "C" int TINKER_MOD(atoms, n);
extern "C" int TINKER_MOD(atoms, type)[maxatm];
extern "C" double TINKER_MOD(atoms, x)[maxatm];
extern "C" double TINKER_MOD(atoms, y)[maxatm];
extern "C" double TINKER_MOD(atoms, z)[maxatm];

int& n = TINKER_MOD(atoms, n);
int (&type)[maxatm] = TINKER_MOD(atoms, type);
double (&x)[maxatm] = TINKER_MOD(atoms, x);
double (&y)[maxatm] = TINKER_MOD(atoms, y);
double (&z)[maxatm] = TINKER_MOD(atoms, z);
#endif

} TINKER_NAMESPACE_END

#endif
