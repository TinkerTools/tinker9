#ifndef TINKER_RT_H_
#define TINKER_RT_H_

#include "cxx.h"

TINKER_NAMESPACE_BEGIN
typedef int logical;
const logical _true_ = 1;
const logical _false_ = 0;

void nextarg(size_t len, char* str, logical& exist);

template <size_t Len>
void nextarg(char (&str)[Len], logical& exist) {
  nextarg(Len, str, exist);
}
TINKER_NAMESPACE_END

extern "C" {
void TINKER_RT(final)();
void TINKER_RT(getxyz)();
void TINKER_RT(initial)();
void TINKER_RT(command)();
void TINKER_RT(mdinit)();
void TINKER_RT(mechanic)();
void TINKER_RT(prterr)();
void TINKER_RT(lattice)();
void TINKER_RT(invert)(int* n, double* a);
void TINKER_RT(mdsave)(int* istep, double* dt, double* epot, double* eksum);

// pmestuf.f
void TINKER_RT(bspline)(double* x, int* n, double* c);
void TINKER_RT(dftmod)(double* bsmod, double* bsarray, int* nfft, int* order);
}

#endif
