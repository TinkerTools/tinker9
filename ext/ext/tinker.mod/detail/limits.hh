#ifndef TINKER_MOD_LIMITS_HH_
#define TINKER_MOD_LIMITS_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace limits {
extern double& vdwcut;
extern double& repcut;
extern double& dispcut;
extern double& chgcut;
extern double& dplcut;
extern double& mpolecut;
extern double& ctrncut;
extern double& vdwtaper;
extern double& reptaper;
extern double& disptaper;
extern double& chgtaper;
extern double& dpltaper;
extern double& mpoletaper;
extern double& ctrntaper;
extern double& ewaldcut;
extern double& dewaldcut;
extern double& usolvcut;
extern int& use_ewald;
extern int& use_dewald;
extern int& use_lights;
extern int& use_list;
extern int& use_vlist;
extern int& use_dlist;
extern int& use_clist;
extern int& use_mlist;
extern int& use_ulist;

#ifdef TINKER_MOD_CPP_
extern "C" double m_tinker_mod(limits, vdwcut);
extern "C" double m_tinker_mod(limits, repcut);
extern "C" double m_tinker_mod(limits, dispcut);
extern "C" double m_tinker_mod(limits, chgcut);
extern "C" double m_tinker_mod(limits, dplcut);
extern "C" double m_tinker_mod(limits, mpolecut);
extern "C" double m_tinker_mod(limits, ctrncut);
extern "C" double m_tinker_mod(limits, vdwtaper);
extern "C" double m_tinker_mod(limits, reptaper);
extern "C" double m_tinker_mod(limits, disptaper);
extern "C" double m_tinker_mod(limits, chgtaper);
extern "C" double m_tinker_mod(limits, dpltaper);
extern "C" double m_tinker_mod(limits, mpoletaper);
extern "C" double m_tinker_mod(limits, ctrntaper);
extern "C" double m_tinker_mod(limits, ewaldcut);
extern "C" double m_tinker_mod(limits, dewaldcut);
extern "C" double m_tinker_mod(limits, usolvcut);
extern "C" int m_tinker_mod(limits, use_ewald);
extern "C" int m_tinker_mod(limits, use_dewald);
extern "C" int m_tinker_mod(limits, use_lights);
extern "C" int m_tinker_mod(limits, use_list);
extern "C" int m_tinker_mod(limits, use_vlist);
extern "C" int m_tinker_mod(limits, use_dlist);
extern "C" int m_tinker_mod(limits, use_clist);
extern "C" int m_tinker_mod(limits, use_mlist);
extern "C" int m_tinker_mod(limits, use_ulist);

double& vdwcut = m_tinker_mod(limits, vdwcut);
double& repcut = m_tinker_mod(limits, repcut);
double& dispcut = m_tinker_mod(limits, dispcut);
double& chgcut = m_tinker_mod(limits, chgcut);
double& dplcut = m_tinker_mod(limits, dplcut);
double& mpolecut = m_tinker_mod(limits, mpolecut);
double& ctrncut = m_tinker_mod(limits, ctrncut);
double& vdwtaper = m_tinker_mod(limits, vdwtaper);
double& reptaper = m_tinker_mod(limits, reptaper);
double& disptaper = m_tinker_mod(limits, disptaper);
double& chgtaper = m_tinker_mod(limits, chgtaper);
double& dpltaper = m_tinker_mod(limits, dpltaper);
double& mpoletaper = m_tinker_mod(limits, mpoletaper);
double& ctrntaper = m_tinker_mod(limits, ctrntaper);
double& ewaldcut = m_tinker_mod(limits, ewaldcut);
double& dewaldcut = m_tinker_mod(limits, dewaldcut);
double& usolvcut = m_tinker_mod(limits, usolvcut);
int& use_ewald = m_tinker_mod(limits, use_ewald);
int& use_dewald = m_tinker_mod(limits, use_dewald);
int& use_lights = m_tinker_mod(limits, use_lights);
int& use_list = m_tinker_mod(limits, use_list);
int& use_vlist = m_tinker_mod(limits, use_vlist);
int& use_dlist = m_tinker_mod(limits, use_dlist);
int& use_clist = m_tinker_mod(limits, use_clist);
int& use_mlist = m_tinker_mod(limits, use_mlist);
int& use_ulist = m_tinker_mod(limits, use_ulist);
#endif

} TINKER_NAMESPACE_END

#endif
