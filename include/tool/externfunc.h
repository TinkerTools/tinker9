#pragma once
#include "tool/error.h"
#include "tool/macro.h"

#define TINKER_F1EXTN(RETURN_TYPE, FUNC, SUFFIX, HAS_IMPL, ...)                                    \
   TINKER_F1EXTN_##HAS_IMPL##__(RETURN_TYPE, FUNC, SUFFIX, __VA_ARGS__)
#define TINKER_F1EXTN_1__(R, F, S, ...) extern R F##_##S(__VA_ARGS__)
#define TINKER_F1EXTN_0__(R, F, S, ...)                                                            \
   static R F##_##S(__VA_ARGS__)                                                                   \
   {                                                                                               \
      throw FatalError("F1EXTN: " #R " " #F "_" #S "(...) is not implemented.");                   \
   }

#define TINKER_F1CALL(FUNC, SUFFIX1, ...) return FUNC##_##SUFFIX1(__VA_ARGS__)

#define TINKER_F2EXTN(RETURN_TYPE, FUNC, SUFFIX1, HAS_IMPL1, SUFFIX2, HAS_IMPL2, ...)              \
   TINKER_F1EXTN(RETURN_TYPE, FUNC, SUFFIX1, HAS_IMPL1, __VA_ARGS__);                              \
   TINKER_F1EXTN(RETURN_TYPE, FUNC, SUFFIX2, HAS_IMPL2, __VA_ARGS__)

#define TINKER_F2CALL(FUNC, SUFFIX1, SUFFIX2, ...)                                                 \
   if (0) {                                                                                        \
   } else if (TINKER_F2CALL_##SUFFIX1##__) {                                                       \
      return FUNC##_##SUFFIX1(__VA_ARGS__);                                                        \
   } else if (TINKER_F2CALL_##SUFFIX2##__) {                                                       \
      return FUNC##_##SUFFIX2(__VA_ARGS__);                                                        \
   } else if (TINKER_HOST) {                                                                       \
      return FUNC##_acc(__VA_ARGS__);                                                              \
   } else {                                                                                        \
      throw FatalError("F2CALL: " #FUNC "(...) does not have a fallback implementation.");         \
   }
#define TINKER_F2CALL_cu__  TINKER_GPULANG_CUDA
#define TINKER_F2CALL_acc__ TINKER_GPULANG_OPENACC
