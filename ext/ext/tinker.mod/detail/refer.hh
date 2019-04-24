#ifndef TINKER_MOD_REFER_HH_
#define TINKER_MOD_REFER_HH_

#include "util/macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace refer {
using namespace sizes;

extern int (&nref)[maxref];
extern int (&refltitle)[maxref];
extern int (&refleng)[maxref];
extern int*& reftyp;
extern int*& n12ref;
extern int*& i12ref;
extern double (&xboxref)[maxref];
extern double (&yboxref)[maxref];
extern double (&zboxref)[maxref];
extern double (&alpharef)[maxref];
extern double (&betaref)[maxref];
extern double (&gammaref)[maxref];
extern double*& xref;
extern double*& yref;
extern double*& zref;
extern char (*&refnam)[3];
extern char (&reffile)[maxref][240];
extern char (&reftitle)[maxref][240];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(refer, nref)[maxref];
extern "C" int m_tinker_mod(refer, refltitle)[maxref];
extern "C" int m_tinker_mod(refer, refleng)[maxref];
extern "C" int* m_tinker_mod(refer, reftyp);
extern "C" int* m_tinker_mod(refer, n12ref);
extern "C" int* m_tinker_mod(refer, i12ref);
extern "C" double m_tinker_mod(refer, xboxref)[maxref];
extern "C" double m_tinker_mod(refer, yboxref)[maxref];
extern "C" double m_tinker_mod(refer, zboxref)[maxref];
extern "C" double m_tinker_mod(refer, alpharef)[maxref];
extern "C" double m_tinker_mod(refer, betaref)[maxref];
extern "C" double m_tinker_mod(refer, gammaref)[maxref];
extern "C" double* m_tinker_mod(refer, xref);
extern "C" double* m_tinker_mod(refer, yref);
extern "C" double* m_tinker_mod(refer, zref);
extern "C" char (*m_tinker_mod(refer, refnam))[3];
extern "C" char m_tinker_mod(refer, reffile)[maxref][240];
extern "C" char m_tinker_mod(refer, reftitle)[maxref][240];

int (&nref)[maxref] = m_tinker_mod(refer, nref);
int (&refltitle)[maxref] = m_tinker_mod(refer, refltitle);
int (&refleng)[maxref] = m_tinker_mod(refer, refleng);
int*& reftyp = m_tinker_mod(refer, reftyp);
int*& n12ref = m_tinker_mod(refer, n12ref);
int*& i12ref = m_tinker_mod(refer, i12ref);
double (&xboxref)[maxref] = m_tinker_mod(refer, xboxref);
double (&yboxref)[maxref] = m_tinker_mod(refer, yboxref);
double (&zboxref)[maxref] = m_tinker_mod(refer, zboxref);
double (&alpharef)[maxref] = m_tinker_mod(refer, alpharef);
double (&betaref)[maxref] = m_tinker_mod(refer, betaref);
double (&gammaref)[maxref] = m_tinker_mod(refer, gammaref);
double*& xref = m_tinker_mod(refer, xref);
double*& yref = m_tinker_mod(refer, yref);
double*& zref = m_tinker_mod(refer, zref);
char (*&refnam)[3] = m_tinker_mod(refer, refnam);
char (&reffile)[maxref][240] = m_tinker_mod(refer, reffile);
char (&reftitle)[maxref][240] = m_tinker_mod(refer, reftitle);
#endif

} TINKER_NAMESPACE_END

#endif
