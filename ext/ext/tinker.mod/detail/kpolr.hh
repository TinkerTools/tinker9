#ifndef TINKER_MOD_KPOLR_HH_
#define TINKER_MOD_KPOLR_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace kpolr {
extern int*& pgrp;
extern double*& polr;
extern double*& athl;

#ifdef TINKER_MOD_CPP_
extern "C" int* m_tinker_mod(kpolr, pgrp);
extern "C" double* m_tinker_mod(kpolr, polr);
extern "C" double* m_tinker_mod(kpolr, athl);

int*& pgrp = m_tinker_mod(kpolr, pgrp);
double*& polr = m_tinker_mod(kpolr, polr);
double*& athl = m_tinker_mod(kpolr, athl);
#endif

} TINKER_NAMESPACE_END

#endif
