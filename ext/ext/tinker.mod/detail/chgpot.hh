#ifndef TINKER_MOD_CHGPOT_HH_
#define TINKER_MOD_CHGPOT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace chgpot {
extern double& electric;
extern double& dielec;
extern double& ebuffer;
extern double& c2scale;
extern double& c3scale;
extern double& c4scale;
extern double& c5scale;
extern int& neutnbr;
extern int& neutcut;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(chgpot, electric);
extern "C" double m_tinker_mod(chgpot, dielec);
extern "C" double m_tinker_mod(chgpot, ebuffer);
extern "C" double m_tinker_mod(chgpot, c2scale);
extern "C" double m_tinker_mod(chgpot, c3scale);
extern "C" double m_tinker_mod(chgpot, c4scale);
extern "C" double m_tinker_mod(chgpot, c5scale);
extern "C" int m_tinker_mod(chgpot, neutnbr);
extern "C" int m_tinker_mod(chgpot, neutcut);

double& electric = m_tinker_mod(chgpot, electric);
double& dielec = m_tinker_mod(chgpot, dielec);
double& ebuffer = m_tinker_mod(chgpot, ebuffer);
double& c2scale = m_tinker_mod(chgpot, c2scale);
double& c3scale = m_tinker_mod(chgpot, c3scale);
double& c4scale = m_tinker_mod(chgpot, c4scale);
double& c5scale = m_tinker_mod(chgpot, c5scale);
int& neutnbr = m_tinker_mod(chgpot, neutnbr);
int& neutcut = m_tinker_mod(chgpot, neutcut);
#endif

} TINKER_NAMESPACE_END

#endif
