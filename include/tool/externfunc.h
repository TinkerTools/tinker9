#pragma once
#include "tool/error.h"
#include "tool/macro.h"

#define TINKER_F1EXTN__(SUFFIX, HAS_IMPL, RETURN_TYPE, FUNC, ...)                                  \
   TINKER_F1EXTN_##HAS_IMPL##__(RETURN_TYPE, FUNC, SUFFIX, __VA_ARGS__)
#define TINKER_F1EXTN_1__(R, F, S, ...) extern R F##_##S(__VA_ARGS__)
#define TINKER_F1EXTN_0__(R, F, S, ...)

#define TINKER_F1CALL(SUFFIX1, FUNC, ...) FUNC##_##SUFFIX1(__VA_ARGS__)

#define TINKER_F2EXTN(SUFFIX1, HAS_IMPL1, SUFFIX2, HAS_IMPL2, RETURN_TYPE, FUNC, ...)              \
   TINKER_F1EXTN__(SUFFIX1, HAS_IMPL1, RETURN_TYPE, FUNC, __VA_ARGS__);                            \
   TINKER_F1EXTN__(SUFFIX2, HAS_IMPL2, RETURN_TYPE, FUNC, __VA_ARGS__)
