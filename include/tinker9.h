#pragma once
#include "macro.h"
#include <string>
#include <tinker/routines.h>

namespace tinker {
void x_analyze(int, char**);
void x_bar(int, char**);
void x_dynamic(int, char**);
void x_help(int, char**);
void x_info(int, char**);
void x_minimize(int, char**);
void x_testgrad(int, char**);
}

#define TINKER9_VERSION_MAJOR 1
#define TINKER9_VERSION_MINOR 0
#define TINKER9_VERSION_PATCH 0

// clang-format off
#ifdef TINKER9_GIT_SHORT_HASH
#   define TINKER9_PROMO1__ "\n" " Commit:       " TINKER_STR(TINKER9_GIT_SHORT_HASH)
#else
#   define TINKER9_PROMO1__ ""
#endif
#ifdef TINKER9_GIT_DATE
#   define TINKER9_PROMO2__ "\n" " Commit Date:  " TINKER9_GIT_DATE
#else
#   define TINKER9_PROMO2__ "\n" " No GIT History"
#endif
#define TINKER9_PROMO3__ " Compiled at:  " __TIME__ "  " __DATE__
#define TINKER9_PROMO_STRING                                                          \
                                                                                 "\n" \
"     ######################################################################    ""\n" \
"   ##########################################################################  ""\n" \
"  ###                                                                      ### ""\n" \
" ###            Tinker9  --  Software Tools for Molecular Design            ###""\n" \
" ##                                                                          ##""\n" \
" ##                      Version 1.0.0-rc  January 2021                      ##""\n" \
" ##                                                                          ##""\n" \
" ##                 Copyright (c)  Zhi Wang & the Ponder Lab                 ##""\n" \
" ###                           All Rights Reserved                          ###""\n" \
"  ###                                                                      ### ""\n" \
"   ##########################################################################  ""\n" \
"     ######################################################################    ""\n" \
                                                                                 "\n" \
TINKER9_PROMO3__                                                                      \
TINKER9_PROMO2__                                                                      \
TINKER9_PROMO1__                                                                 "\n"
// clang-format on

namespace tinker {
void promo();
void initial();
void mechanic2();

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
