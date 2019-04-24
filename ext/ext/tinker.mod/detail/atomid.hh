#ifndef TINKER_MOD_ATOMID_HH_
#define TINKER_MOD_ATOMID_HH_

#include "util/macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace atomid {
using namespace sizes;

extern int (&tag)[maxatm];
extern int (&class_)[maxatm];
extern int (&atomic)[maxatm];
extern int (&valence)[maxatm];
extern double (&mass)[maxatm];
extern char (&name)[maxatm][3];
extern char (&story)[maxatm][24];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(atomid, tag)[maxatm];
extern "C" int m_tinker_mod(atomid, class)[maxatm];
extern "C" int m_tinker_mod(atomid, atomic)[maxatm];
extern "C" int m_tinker_mod(atomid, valence)[maxatm];
extern "C" double m_tinker_mod(atomid, mass)[maxatm];
extern "C" char m_tinker_mod(atomid, name)[maxatm][3];
extern "C" char m_tinker_mod(atomid, story)[maxatm][24];

int (&tag)[maxatm] = m_tinker_mod(atomid, tag);
int (&class_)[maxatm] = m_tinker_mod(atomid, class);
int (&atomic)[maxatm] = m_tinker_mod(atomid, atomic);
int (&valence)[maxatm] = m_tinker_mod(atomid, valence);
double (&mass)[maxatm] = m_tinker_mod(atomid, mass);
char (&name)[maxatm][3] = m_tinker_mod(atomid, name);
char (&story)[maxatm][24] = m_tinker_mod(atomid, story);
#endif

} TINKER_NAMESPACE_END

#endif
