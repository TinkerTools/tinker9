#ifndef TINKER_MOD_OPENMP_HH_
#define TINKER_MOD_OPENMP_HH_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN namespace openmp {
extern int& nproc;
extern int& nthread;

#ifdef TINKER_MOD_CPP_
extern "C" int m_tinker_mod(openmp, nproc);
extern "C" int m_tinker_mod(openmp, nthread);

int& nproc = m_tinker_mod(openmp, nproc);
int& nthread = m_tinker_mod(openmp, nthread);
#endif

} TINKER_NAMESPACE_END

#endif
