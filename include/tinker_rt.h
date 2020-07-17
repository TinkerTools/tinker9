#pragma once
#include "mathfunc.h"
#include "md.h"
#include "subroutine.h"
#include "tool/io_fort_str.h"
#include "tool/io_print.h"
#include "tool/io_read.h"
#include "tool/io_text.h"
#include <tinker/detail/keys.hh>


namespace tinker {
void nextarg(size_t len, char* str, int& exist);


template <size_t Len>
void nextarg(char (&str)[Len], int& exist)
{
   nextarg(Len, str, exist);
}


template <class T1, class T2>
void get_kv(std::string k, T1& v, T2 vdefault);


template <class T>
void get_kbool(std::string k, T& v, bool v_if_k_not_found);
}


extern "C"
{
   void TINKER_RT(final)();
   void TINKER_RT(getxyz)();
   void TINKER_RT(command)();
   void TINKER_RT(mdinit)();
   void TINKER_RT(invert)(int* n, double* a);
   void TINKER_RT(mdsave)(int* istep, double* dt, double* epot, double* eksum);

   // pmestuf.f
   void TINKER_RT(bspline)(double* x, int* n, double* c);
   void TINKER_RT(dftmod)(double* bsmod, double* bsarray, int* nfft,
                          int* order);
}
