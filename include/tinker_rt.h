#pragma once
#include "io_fort_str.h"
#include "io_print.h"
#include "io_read.h"
#include "io_text.h"
#include "mathfunc.h"
#include "md.h"
#include <tinker/detail/keys.hh>


TINKER_NAMESPACE_BEGIN
void nextarg(size_t len, char* str, int& exist);


template <size_t Len>
void nextarg(char (&str)[Len], int& exist)
{
   nextarg(Len, str, exist);
}


template <class T1, class T2>
void get_kv_pair(std::string k, T1& v, T2 vdefault);
TINKER_NAMESPACE_END


extern "C"
{
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
   void TINKER_RT(image)(double*, double*, double*);
   void TINKER_RT(imagen)(double*, double*, double*);

   // pmestuf.f
   void TINKER_RT(bspline)(double* x, int* n, double* c);
   void TINKER_RT(dftmod)(double* bsmod, double* bsarray, int* nfft,
                          int* order);
}
