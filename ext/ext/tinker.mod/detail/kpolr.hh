#ifndef TINKER_MOD_KPOLR_HH_
#define TINKER_MOD_KPOLR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kpolr {
extern int*& pgrp;
extern double*& polr;
extern double*& athl;

#ifdef TINKER_MOD_CPP_
extern "C" int* TINKER_MOD(kpolr, pgrp);
extern "C" double* TINKER_MOD(kpolr, polr);
extern "C" double* TINKER_MOD(kpolr, athl);

int*& pgrp = TINKER_MOD(kpolr, pgrp);
double*& polr = TINKER_MOD(kpolr, polr);
double*& athl = TINKER_MOD(kpolr, athl);
#endif

} TINKER_NAMESPACE_END

#endif
