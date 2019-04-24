#ifndef TINKER_MOD_ENERGI_HH_
#define TINKER_MOD_ENERGI_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace energi {
extern double& esum;
extern double& eb;
extern double& ea;
extern double& eba;
extern double& eub;
extern double& eaa;
extern double& eopb;
extern double& eopd;
extern double& eid;
extern double& eit;
extern double& et;
extern double& ept;
extern double& ebt;
extern double& eat;
extern double& ett;
extern double& ev;
extern double& er;
extern double& edsp;
extern double& ec;
extern double& ecd;
extern double& ed;
extern double& em;
extern double& ep;
extern double& ect;
extern double& erxf;
extern double& es;
extern double& elf;
extern double& eg;
extern double& ex;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(energi, esum);
extern "C" double m_tinker_mod(energi, eb);
extern "C" double m_tinker_mod(energi, ea);
extern "C" double m_tinker_mod(energi, eba);
extern "C" double m_tinker_mod(energi, eub);
extern "C" double m_tinker_mod(energi, eaa);
extern "C" double m_tinker_mod(energi, eopb);
extern "C" double m_tinker_mod(energi, eopd);
extern "C" double m_tinker_mod(energi, eid);
extern "C" double m_tinker_mod(energi, eit);
extern "C" double m_tinker_mod(energi, et);
extern "C" double m_tinker_mod(energi, ept);
extern "C" double m_tinker_mod(energi, ebt);
extern "C" double m_tinker_mod(energi, eat);
extern "C" double m_tinker_mod(energi, ett);
extern "C" double m_tinker_mod(energi, ev);
extern "C" double m_tinker_mod(energi, er);
extern "C" double m_tinker_mod(energi, edsp);
extern "C" double m_tinker_mod(energi, ec);
extern "C" double m_tinker_mod(energi, ecd);
extern "C" double m_tinker_mod(energi, ed);
extern "C" double m_tinker_mod(energi, em);
extern "C" double m_tinker_mod(energi, ep);
extern "C" double m_tinker_mod(energi, ect);
extern "C" double m_tinker_mod(energi, erxf);
extern "C" double m_tinker_mod(energi, es);
extern "C" double m_tinker_mod(energi, elf);
extern "C" double m_tinker_mod(energi, eg);
extern "C" double m_tinker_mod(energi, ex);

double& esum = m_tinker_mod(energi, esum);
double& eb = m_tinker_mod(energi, eb);
double& ea = m_tinker_mod(energi, ea);
double& eba = m_tinker_mod(energi, eba);
double& eub = m_tinker_mod(energi, eub);
double& eaa = m_tinker_mod(energi, eaa);
double& eopb = m_tinker_mod(energi, eopb);
double& eopd = m_tinker_mod(energi, eopd);
double& eid = m_tinker_mod(energi, eid);
double& eit = m_tinker_mod(energi, eit);
double& et = m_tinker_mod(energi, et);
double& ept = m_tinker_mod(energi, ept);
double& ebt = m_tinker_mod(energi, ebt);
double& eat = m_tinker_mod(energi, eat);
double& ett = m_tinker_mod(energi, ett);
double& ev = m_tinker_mod(energi, ev);
double& er = m_tinker_mod(energi, er);
double& edsp = m_tinker_mod(energi, edsp);
double& ec = m_tinker_mod(energi, ec);
double& ecd = m_tinker_mod(energi, ecd);
double& ed = m_tinker_mod(energi, ed);
double& em = m_tinker_mod(energi, em);
double& ep = m_tinker_mod(energi, ep);
double& ect = m_tinker_mod(energi, ect);
double& erxf = m_tinker_mod(energi, erxf);
double& es = m_tinker_mod(energi, es);
double& elf = m_tinker_mod(energi, elf);
double& eg = m_tinker_mod(energi, eg);
double& ex = m_tinker_mod(energi, ex);
#endif

} TINKER_NAMESPACE_END

#endif
