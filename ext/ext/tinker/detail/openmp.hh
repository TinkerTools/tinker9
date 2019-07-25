#ifndef TINKER_MOD_OPENMP_HH_
#define TINKER_MOD_OPENMP_HH_

#include "util_macro.h"

TINKER_NAMESPACE_BEGIN namespace openmp {
extern int& nproc;
extern int& nthread;

#ifdef TINKER_MOD_CPP_
extern "C" int TINKER_MOD(openmp, nproc);
extern "C" int TINKER_MOD(openmp, nthread);

int& nproc = TINKER_MOD(openmp, nproc);
int& nthread = TINKER_MOD(openmp, nthread);
#endif
} TINKER_NAMESPACE_END

#endif
