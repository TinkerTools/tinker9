#pragma once

#include "macro.h"

TINKER_NAMESPACE_BEGIN namespace kpolr {
extern int*& pgrp;
extern double*& polr;
extern double*& athl;
extern double*& ddir;

#ifdef TINKER_MOD_CPP_
extern "C" int* TINKER_MOD(kpolr, pgrp);
extern "C" double* TINKER_MOD(kpolr, polr);
extern "C" double* TINKER_MOD(kpolr, athl);
extern "C" double* TINKER_MOD(kpolr, ddir);

int*& pgrp = TINKER_MOD(kpolr, pgrp);
double*& polr = TINKER_MOD(kpolr, polr);
double*& athl = TINKER_MOD(kpolr, athl);
double*& ddir = TINKER_MOD(kpolr, ddir);
#endif
} TINKER_NAMESPACE_END
