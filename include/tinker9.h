#pragma once
#include "tool/macro.h"
#include <string>

/// \def TINKER9_VERSION_MAJOR
/// \ingroup general
/// \def TINKER9_VERSION_MINOR
/// \ingroup general
/// \def TINKER9_VERSION_PATCH
/// \ingroup general
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
" ##                       Version 1.0.0-rc  April 2021                       ##""\n" \
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
/// \ingroup general
/// \{
void xAnalyze(int, char**);  ///< Entry point of the \c analyze program.
void xBar(int, char**);      ///< Entry point of the \c bar program.
void xDynamic(int, char**);  ///< Entry point of the \c dynamic program.
void xHelp(int, char**);     ///< Entry point of the \c help program.
void xInfo(int, char**);     ///< Entry point of the \c info program.
void xMinimize(int, char**); ///< Entry point of the \c minimize program.
void xTestgrad(int, char**); ///< Entry point of the \c testgrad program.
/// \}
}

namespace tinker {
/// \ingroup general
/// \{
/// \brief Writes a banner message.
void promo();
/// \brief Sets up original values. This function must be translated
/// from the Fortran \c initial subroutine line-by-line.
void initial();
/// \brief Sets up extra parameters and selectable options
/// in addition to the Fortran \c mechanic subroutine.
void mechanic2();
/// \}
}
