#ifndef TINKER_MOD_PTABLE_HH_
#define TINKER_MOD_PTABLE_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace ptable {
const int maxele = 112;
extern double (&atmass)[maxele];
extern double (&vdwrad)[maxele];
extern double (&covrad)[maxele];
extern char (&elemnt)[maxele][3];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(ptable, atmass)[maxele];
extern "C" double m_tinker_mod(ptable, vdwrad)[maxele];
extern "C" double m_tinker_mod(ptable, covrad)[maxele];
extern "C" char m_tinker_mod(ptable, elemnt)[maxele][3];

double (&atmass)[maxele] = m_tinker_mod(ptable, atmass);
double (&vdwrad)[maxele] = m_tinker_mod(ptable, vdwrad);
double (&covrad)[maxele] = m_tinker_mod(ptable, covrad);
char (&elemnt)[maxele][3] = m_tinker_mod(ptable, elemnt);
#endif

} TINKER_NAMESPACE_END

#endif
