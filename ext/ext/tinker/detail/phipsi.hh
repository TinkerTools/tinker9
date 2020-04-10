#ifndef TINKER_MOD_PHIPSI_HH_
#define TINKER_MOD_PHIPSI_HH_

#include "macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace phipsi {
using namespace sizes;

extern int (&chiral)[maxres];
extern int (&disulf)[maxres];
extern double (&phi)[maxres];
extern double (&psi)[maxres];
extern double (&omg)[maxres];
extern double (&chi)[maxres][4];

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(phipsi, chiral)[maxres];
extern "C" int TINKER_MOD(phipsi, disulf)[maxres];
extern "C" double TINKER_MOD(phipsi, phi)[maxres];
extern "C" double TINKER_MOD(phipsi, psi)[maxres];
extern "C" double TINKER_MOD(phipsi, omg)[maxres];
extern "C" double TINKER_MOD(phipsi, chi)[maxres][4];

int (&chiral)[maxres] = TINKER_MOD(phipsi, chiral);
int (&disulf)[maxres] = TINKER_MOD(phipsi, disulf);
double (&phi)[maxres] = TINKER_MOD(phipsi, phi);
double (&psi)[maxres] = TINKER_MOD(phipsi, psi);
double (&omg)[maxres] = TINKER_MOD(phipsi, omg);
double (&chi)[maxres][4] = TINKER_MOD(phipsi, chi);
#endif
} TINKER_NAMESPACE_END

#endif
