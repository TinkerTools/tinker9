#ifndef TINKER_MOD_PHIPSI_HH_
#define TINKER_MOD_PHIPSI_HH_

#include "util/macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace phipsi {
using namespace sizes;

extern int (&chiral)[maxres];
extern int (&disulf)[maxres];
extern double (&phi)[maxres];
extern double (&psi)[maxres];
extern double (&omega)[maxres];
extern double (&chi)[maxres][4];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(phipsi, chiral)[maxres];
extern "C" int m_tinker_mod(phipsi, disulf)[maxres];
extern "C" double m_tinker_mod(phipsi, phi)[maxres];
extern "C" double m_tinker_mod(phipsi, psi)[maxres];
extern "C" double m_tinker_mod(phipsi, omega)[maxres];
extern "C" double m_tinker_mod(phipsi, chi)[maxres][4];

int (&chiral)[maxres] = m_tinker_mod(phipsi, chiral);
int (&disulf)[maxres] = m_tinker_mod(phipsi, disulf);
double (&phi)[maxres] = m_tinker_mod(phipsi, phi);
double (&psi)[maxres] = m_tinker_mod(phipsi, psi);
double (&omega)[maxres] = m_tinker_mod(phipsi, omega);
double (&chi)[maxres][4] = m_tinker_mod(phipsi, chi);
#endif

} TINKER_NAMESPACE_END

#endif
