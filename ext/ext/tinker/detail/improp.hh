#ifndef TINKER_MOD_IMPROP_HH_
#define TINKER_MOD_IMPROP_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace improp {
extern int& niprop;
extern int*& iiprop;
extern double*& kprop;
extern double*& vprop;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(improp, niprop);
extern "C" int* TINKER_MOD(improp, iiprop);
extern "C" double* TINKER_MOD(improp, kprop);
extern "C" double* TINKER_MOD(improp, vprop);

int& niprop = TINKER_MOD(improp, niprop);
int*& iiprop = TINKER_MOD(improp, iiprop);
double*& kprop = TINKER_MOD(improp, kprop);
double*& vprop = TINKER_MOD(improp, vprop);
#endif

} TINKER_NAMESPACE_END

#endif
