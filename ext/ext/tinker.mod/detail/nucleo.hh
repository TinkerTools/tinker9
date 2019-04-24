#ifndef TINKER_MOD_NUCLEO_HH_
#define TINKER_MOD_NUCLEO_HH_

#include "util/macro.h"
#include "sizes.hh"

TINKER_NAMESPACE_BEGIN namespace nucleo {
using namespace sizes;

extern int (&pucker)[maxres];
extern double (&glyco)[maxres];
extern double (&bkbone)[maxres][6];
extern int& dblhlx;
extern int (&deoxy)[maxres];
extern char (&hlxform)[1];

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(nucleo, pucker)[maxres];
extern "C" double m_tinker_mod(nucleo, glyco)[maxres];
extern "C" double m_tinker_mod(nucleo, bkbone)[maxres][6];
extern "C" int m_tinker_mod(nucleo, dblhlx);
extern "C" int m_tinker_mod(nucleo, deoxy)[maxres];
extern "C" char m_tinker_mod(nucleo, hlxform)[1];

int (&pucker)[maxres] = m_tinker_mod(nucleo, pucker);
double (&glyco)[maxres] = m_tinker_mod(nucleo, glyco);
double (&bkbone)[maxres][6] = m_tinker_mod(nucleo, bkbone);
int& dblhlx = m_tinker_mod(nucleo, dblhlx);
int (&deoxy)[maxres] = m_tinker_mod(nucleo, deoxy);
char (&hlxform)[1] = m_tinker_mod(nucleo, hlxform);
#endif

} TINKER_NAMESPACE_END

#endif
