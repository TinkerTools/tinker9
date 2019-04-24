#ifndef TINKER_MOD_SHUNT_HH_
#define TINKER_MOD_SHUNT_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace shunt {
extern double& off;
extern double& off2;
extern double& cut;
extern double& cut2;
extern double& c0;
extern double& c1;
extern double& c2;
extern double& c3;
extern double& c4;
extern double& c5;
extern double& f0;
extern double& f1;
extern double& f2;
extern double& f3;
extern double& f4;
extern double& f5;
extern double& f6;
extern double& f7;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(shunt, off);
extern "C" double m_tinker_mod(shunt, off2);
extern "C" double m_tinker_mod(shunt, cut);
extern "C" double m_tinker_mod(shunt, cut2);
extern "C" double m_tinker_mod(shunt, c0);
extern "C" double m_tinker_mod(shunt, c1);
extern "C" double m_tinker_mod(shunt, c2);
extern "C" double m_tinker_mod(shunt, c3);
extern "C" double m_tinker_mod(shunt, c4);
extern "C" double m_tinker_mod(shunt, c5);
extern "C" double m_tinker_mod(shunt, f0);
extern "C" double m_tinker_mod(shunt, f1);
extern "C" double m_tinker_mod(shunt, f2);
extern "C" double m_tinker_mod(shunt, f3);
extern "C" double m_tinker_mod(shunt, f4);
extern "C" double m_tinker_mod(shunt, f5);
extern "C" double m_tinker_mod(shunt, f6);
extern "C" double m_tinker_mod(shunt, f7);

double& off = m_tinker_mod(shunt, off);
double& off2 = m_tinker_mod(shunt, off2);
double& cut = m_tinker_mod(shunt, cut);
double& cut2 = m_tinker_mod(shunt, cut2);
double& c0 = m_tinker_mod(shunt, c0);
double& c1 = m_tinker_mod(shunt, c1);
double& c2 = m_tinker_mod(shunt, c2);
double& c3 = m_tinker_mod(shunt, c3);
double& c4 = m_tinker_mod(shunt, c4);
double& c5 = m_tinker_mod(shunt, c5);
double& f0 = m_tinker_mod(shunt, f0);
double& f1 = m_tinker_mod(shunt, f1);
double& f2 = m_tinker_mod(shunt, f2);
double& f3 = m_tinker_mod(shunt, f3);
double& f4 = m_tinker_mod(shunt, f4);
double& f5 = m_tinker_mod(shunt, f5);
double& f6 = m_tinker_mod(shunt, f6);
double& f7 = m_tinker_mod(shunt, f7);
#endif

} TINKER_NAMESPACE_END

#endif
