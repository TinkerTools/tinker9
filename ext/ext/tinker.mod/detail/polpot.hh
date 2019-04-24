#ifndef TINKER_MOD_POLPOT_HH_
#define TINKER_MOD_POLPOT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace polpot {
extern int& politer;
extern double& poleps;
extern double& p2scale;
extern double& p3scale;
extern double& p4scale;
extern double& p5scale;
extern double& p2iscale;
extern double& p3iscale;
extern double& p4iscale;
extern double& p5iscale;
extern double& d1scale;
extern double& d2scale;
extern double& d3scale;
extern double& d4scale;
extern double& u1scale;
extern double& u2scale;
extern double& u3scale;
extern double& u4scale;
extern double& w2scale;
extern double& w3scale;
extern double& w4scale;
extern double& w5scale;
extern double& udiag;
extern int& dpequal;
extern int& use_thole;
extern char (&poltyp)[6];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(polpot, politer);
extern "C" double m_tinker_mod(polpot, poleps);
extern "C" double m_tinker_mod(polpot, p2scale);
extern "C" double m_tinker_mod(polpot, p3scale);
extern "C" double m_tinker_mod(polpot, p4scale);
extern "C" double m_tinker_mod(polpot, p5scale);
extern "C" double m_tinker_mod(polpot, p2iscale);
extern "C" double m_tinker_mod(polpot, p3iscale);
extern "C" double m_tinker_mod(polpot, p4iscale);
extern "C" double m_tinker_mod(polpot, p5iscale);
extern "C" double m_tinker_mod(polpot, d1scale);
extern "C" double m_tinker_mod(polpot, d2scale);
extern "C" double m_tinker_mod(polpot, d3scale);
extern "C" double m_tinker_mod(polpot, d4scale);
extern "C" double m_tinker_mod(polpot, u1scale);
extern "C" double m_tinker_mod(polpot, u2scale);
extern "C" double m_tinker_mod(polpot, u3scale);
extern "C" double m_tinker_mod(polpot, u4scale);
extern "C" double m_tinker_mod(polpot, w2scale);
extern "C" double m_tinker_mod(polpot, w3scale);
extern "C" double m_tinker_mod(polpot, w4scale);
extern "C" double m_tinker_mod(polpot, w5scale);
extern "C" double m_tinker_mod(polpot, udiag);
extern "C" int m_tinker_mod(polpot, dpequal);
extern "C" int m_tinker_mod(polpot, use_thole);
extern "C" char m_tinker_mod(polpot, poltyp)[6];

int& politer = m_tinker_mod(polpot, politer);
double& poleps = m_tinker_mod(polpot, poleps);
double& p2scale = m_tinker_mod(polpot, p2scale);
double& p3scale = m_tinker_mod(polpot, p3scale);
double& p4scale = m_tinker_mod(polpot, p4scale);
double& p5scale = m_tinker_mod(polpot, p5scale);
double& p2iscale = m_tinker_mod(polpot, p2iscale);
double& p3iscale = m_tinker_mod(polpot, p3iscale);
double& p4iscale = m_tinker_mod(polpot, p4iscale);
double& p5iscale = m_tinker_mod(polpot, p5iscale);
double& d1scale = m_tinker_mod(polpot, d1scale);
double& d2scale = m_tinker_mod(polpot, d2scale);
double& d3scale = m_tinker_mod(polpot, d3scale);
double& d4scale = m_tinker_mod(polpot, d4scale);
double& u1scale = m_tinker_mod(polpot, u1scale);
double& u2scale = m_tinker_mod(polpot, u2scale);
double& u3scale = m_tinker_mod(polpot, u3scale);
double& u4scale = m_tinker_mod(polpot, u4scale);
double& w2scale = m_tinker_mod(polpot, w2scale);
double& w3scale = m_tinker_mod(polpot, w3scale);
double& w4scale = m_tinker_mod(polpot, w4scale);
double& w5scale = m_tinker_mod(polpot, w5scale);
double& udiag = m_tinker_mod(polpot, udiag);
int& dpequal = m_tinker_mod(polpot, dpequal);
int& use_thole = m_tinker_mod(polpot, use_thole);
char (&poltyp)[6] = m_tinker_mod(polpot, poltyp);
#endif

} TINKER_NAMESPACE_END

#endif
