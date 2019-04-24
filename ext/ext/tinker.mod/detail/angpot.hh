#ifndef TINKER_MOD_ANGPOT_HH_
#define TINKER_MOD_ANGPOT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace angpot {
extern double& angunit;
extern double& stbnunit;
extern double& aaunit;
extern double& opbunit;
extern double& opdunit;
extern double& cang;
extern double& qang;
extern double& pang;
extern double& sang;
extern double& copb;
extern double& qopb;
extern double& popb;
extern double& sopb;
extern double& copd;
extern double& qopd;
extern double& popd;
extern double& sopd;
extern char (&angtrig)[8];
extern char (&opbtyp)[8];
extern char (*&angtyp)[8];

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(angpot, angunit);
extern "C" double m_tinker_mod(angpot, stbnunit);
extern "C" double m_tinker_mod(angpot, aaunit);
extern "C" double m_tinker_mod(angpot, opbunit);
extern "C" double m_tinker_mod(angpot, opdunit);
extern "C" double m_tinker_mod(angpot, cang);
extern "C" double m_tinker_mod(angpot, qang);
extern "C" double m_tinker_mod(angpot, pang);
extern "C" double m_tinker_mod(angpot, sang);
extern "C" double m_tinker_mod(angpot, copb);
extern "C" double m_tinker_mod(angpot, qopb);
extern "C" double m_tinker_mod(angpot, popb);
extern "C" double m_tinker_mod(angpot, sopb);
extern "C" double m_tinker_mod(angpot, copd);
extern "C" double m_tinker_mod(angpot, qopd);
extern "C" double m_tinker_mod(angpot, popd);
extern "C" double m_tinker_mod(angpot, sopd);
extern "C" char m_tinker_mod(angpot, angtrig)[8];
extern "C" char m_tinker_mod(angpot, opbtyp)[8];
extern "C" char (*m_tinker_mod(angpot, angtyp))[8];

double& angunit = m_tinker_mod(angpot, angunit);
double& stbnunit = m_tinker_mod(angpot, stbnunit);
double& aaunit = m_tinker_mod(angpot, aaunit);
double& opbunit = m_tinker_mod(angpot, opbunit);
double& opdunit = m_tinker_mod(angpot, opdunit);
double& cang = m_tinker_mod(angpot, cang);
double& qang = m_tinker_mod(angpot, qang);
double& pang = m_tinker_mod(angpot, pang);
double& sang = m_tinker_mod(angpot, sang);
double& copb = m_tinker_mod(angpot, copb);
double& qopb = m_tinker_mod(angpot, qopb);
double& popb = m_tinker_mod(angpot, popb);
double& sopb = m_tinker_mod(angpot, sopb);
double& copd = m_tinker_mod(angpot, copd);
double& qopd = m_tinker_mod(angpot, qopd);
double& popd = m_tinker_mod(angpot, popd);
double& sopd = m_tinker_mod(angpot, sopd);
char (&angtrig)[8] = m_tinker_mod(angpot, angtrig);
char (&opbtyp)[8] = m_tinker_mod(angpot, opbtyp);
char (*&angtyp)[8] = m_tinker_mod(angpot, angtyp);
#endif

} TINKER_NAMESPACE_END

#endif
