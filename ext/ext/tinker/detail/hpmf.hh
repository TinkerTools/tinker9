#ifndef TINKER_MOD_HPMF_HH_
#define TINKER_MOD_HPMF_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace hpmf {
const double rcarbon = 1.7e0;
const double rwater = 1.4e0;
const double acsurf = 120.7628e0;
const double safact = 0.3516e0;
const double tslope = 100.0e0;
const double toffset = 6.0e0;
const double hpmfcut = 11.0e0;
const double h1 = -0.7308004860404441194e0;
const double h2 = 0.2001645051578760659e0;
const double h3 = -0.0905499953418473502e0;
const double c1 = 3.8167879266271396155e0;
const double c2 = 5.4669162286016419472e0;
const double c3 = 7.1167694861385353278e0;
const double w1 = 1.6858993102248638341e0;
const double w2 = 1.3906405621629980285e0;
const double w3 = 1.5741657341338335385e0;
extern int& npmf;
extern int*& ipmf;
extern double*& rpmf;
extern double*& acsa;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(hpmf, npmf);
extern "C" int* TINKER_MOD(hpmf, ipmf);
extern "C" double* TINKER_MOD(hpmf, rpmf);
extern "C" double* TINKER_MOD(hpmf, acsa);

int& npmf = TINKER_MOD(hpmf, npmf);
int*& ipmf = TINKER_MOD(hpmf, ipmf);
double*& rpmf = TINKER_MOD(hpmf, rpmf);
double*& acsa = TINKER_MOD(hpmf, acsa);
#endif
} TINKER_NAMESPACE_END

#endif
