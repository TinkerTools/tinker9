#ifndef TINKER_MOD_NUCLEO_HH_
#define TINKER_MOD_NUCLEO_HH_

#include "macro.h"
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
extern "C" int TINKER_MOD(nucleo, pucker)[maxres];
extern "C" double TINKER_MOD(nucleo, glyco)[maxres];
extern "C" double TINKER_MOD(nucleo, bkbone)[maxres][6];
extern "C" int TINKER_MOD(nucleo, dblhlx);
extern "C" int TINKER_MOD(nucleo, deoxy)[maxres];
extern "C" char TINKER_MOD(nucleo, hlxform)[1];

int (&pucker)[maxres] = TINKER_MOD(nucleo, pucker);
double (&glyco)[maxres] = TINKER_MOD(nucleo, glyco);
double (&bkbone)[maxres][6] = TINKER_MOD(nucleo, bkbone);
int& dblhlx = TINKER_MOD(nucleo, dblhlx);
int (&deoxy)[maxres] = TINKER_MOD(nucleo, deoxy);
char (&hlxform)[1] = TINKER_MOD(nucleo, hlxform);
#endif
} TINKER_NAMESPACE_END

#endif
