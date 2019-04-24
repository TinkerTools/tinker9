#ifndef TINKER_MOD_VDWPOT_HH_
#define TINKER_MOD_VDWPOT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace vdwpot {
const int maxgauss = 10;
extern int& ngauss;
extern double (&igauss)[maxgauss][2];
extern double& abuck;
extern double& bbuck;
extern double& cbuck;
extern double& ghal;
extern double& dhal;
extern double& v2scale;
extern double& v3scale;
extern double& v4scale;
extern double& v5scale;
extern int& use_vcorr;
extern char (&vdwindex)[5];
extern char (&radtyp)[5];
extern char (&radsiz)[8];
extern char (&gausstyp)[8];
extern char (&radrule)[10];
extern char (&epsrule)[10];
extern char (&vdwtyp)[13];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(vdwpot, ngauss);
extern "C" double m_tinker_mod(vdwpot, igauss)[maxgauss][2];
extern "C" double m_tinker_mod(vdwpot, abuck);
extern "C" double m_tinker_mod(vdwpot, bbuck);
extern "C" double m_tinker_mod(vdwpot, cbuck);
extern "C" double m_tinker_mod(vdwpot, ghal);
extern "C" double m_tinker_mod(vdwpot, dhal);
extern "C" double m_tinker_mod(vdwpot, v2scale);
extern "C" double m_tinker_mod(vdwpot, v3scale);
extern "C" double m_tinker_mod(vdwpot, v4scale);
extern "C" double m_tinker_mod(vdwpot, v5scale);
extern "C" int m_tinker_mod(vdwpot, use_vcorr);
extern "C" char m_tinker_mod(vdwpot, vdwindex)[5];
extern "C" char m_tinker_mod(vdwpot, radtyp)[5];
extern "C" char m_tinker_mod(vdwpot, radsiz)[8];
extern "C" char m_tinker_mod(vdwpot, gausstyp)[8];
extern "C" char m_tinker_mod(vdwpot, radrule)[10];
extern "C" char m_tinker_mod(vdwpot, epsrule)[10];
extern "C" char m_tinker_mod(vdwpot, vdwtyp)[13];

int& ngauss = m_tinker_mod(vdwpot, ngauss);
double (&igauss)[maxgauss][2] = m_tinker_mod(vdwpot, igauss);
double& abuck = m_tinker_mod(vdwpot, abuck);
double& bbuck = m_tinker_mod(vdwpot, bbuck);
double& cbuck = m_tinker_mod(vdwpot, cbuck);
double& ghal = m_tinker_mod(vdwpot, ghal);
double& dhal = m_tinker_mod(vdwpot, dhal);
double& v2scale = m_tinker_mod(vdwpot, v2scale);
double& v3scale = m_tinker_mod(vdwpot, v3scale);
double& v4scale = m_tinker_mod(vdwpot, v4scale);
double& v5scale = m_tinker_mod(vdwpot, v5scale);
int& use_vcorr = m_tinker_mod(vdwpot, use_vcorr);
char (&vdwindex)[5] = m_tinker_mod(vdwpot, vdwindex);
char (&radtyp)[5] = m_tinker_mod(vdwpot, radtyp);
char (&radsiz)[8] = m_tinker_mod(vdwpot, radsiz);
char (&gausstyp)[8] = m_tinker_mod(vdwpot, gausstyp);
char (&radrule)[10] = m_tinker_mod(vdwpot, radrule);
char (&epsrule)[10] = m_tinker_mod(vdwpot, epsrule);
char (&vdwtyp)[13] = m_tinker_mod(vdwpot, vdwtyp);
#endif

} TINKER_NAMESPACE_END

#endif
