#pragma once

namespace tinker {
void promo();
}

#define TINKER9_VERSION_MAJOR 1
#define TINKER9_VERSION_MINOR 0
#define TINKER9_VERSION_PATCH 0

#ifdef TINKER9_GIT_SHORT_HASH
#   define TINKER9_PROMO1__                                                                        \
      "\n"                                                                                         \
      " Commit:       " TINKER_STR(TINKER9_GIT_SHORT_HASH)
#else
#   define TINKER9_PROMO1__ ""
#endif
#ifdef TINKER9_GIT_DATE
#   define TINKER9_PROMO2__                                                                        \
      "\n"                                                                                         \
      " Commit Date:  " TINKER9_GIT_DATE
#else
#   define TINKER9_PROMO2__                                                                        \
      "\n"                                                                                         \
      " No GIT History"
#endif
#define TINKER9_PROMO3__ " Compiled at:  " __TIME__ "  " __DATE__
// clang-format off
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
