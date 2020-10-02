#pragma once


/**
 * \def TINKER9_VERSION_MAJOR
 * \ingroup macro
 * Major version.
 *
 * \def TINKER9_VERSION_MINOR
 * \ingroup macro
 * Minor version.
 *
 * \def TINKER9_VERSION_PATCH
 * \ingroup macro
 * Patch version.
 */
#define TINKER9_VERSION_MAJOR 0
#define TINKER9_VERSION_MINOR 3
#define TINKER9_VERSION_PATCH 30


#ifdef TINKER9_GIT_SHORT_HASH
#   define TINKER9_PROMO1__ " Commit " TINKER_STR(TINKER9_GIT_SHORT_HASH)
#else
#   define TINKER9_PROMO1__ ""
#endif
#define TINKER9_PROMO2__ " Compiled at " __TIME__ "  " __DATE__
// clang-format off
/**
 * \def TINKER9_PROMO_STRING
 * \ingroup macro
 */
#define TINKER9_PROMO_STRING                                                   \
                                                                           "\n"\
"     ###############################################################    " "\n"\
"   ###################################################################  " "\n"\
"  ###                                                               ### " "\n"\
" ###        Tinker9  ---  Software Tools for Molecular Design        ###" "\n"\
" ##                                                                   ##" "\n"\
" ##                      Alpha Testing  Oct 2020                      ##" "\n"\
" ###                       All Rights Reserved                       ###" "\n"\
"  ###                                                               ### " "\n"\
"   ###################################################################  " "\n"\
"     ###############################################################    " "\n"\
TINKER9_PROMO1__ TINKER9_PROMO2__                                "\n"
// clang-format on
