#ifndef TINKER_RT_H_
#define TINKER_RT_H_

#include "util/cxx.h"

TINKER_NAMESPACE_BEGIN
typedef int logical;
const logical _true_ = 1;
const logical _false_ = 0;

void nextarg(size_t _len, char* _str, logical& _exist);

template <size_t _Len>
void nextarg(char (&_str)[_Len], logical& _exist) {
  nextarg(_Len, _str, _exist);
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

// pmestuf.f
void TINKER_RT(bspline)(double* x, int* n, double* c);
void TINKER_RT(dftmod)(double* bsmod, double* bsarray, int* nfft, int* order);
}

#endif
