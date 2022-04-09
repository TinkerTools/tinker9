#pragma once
#include "tool/error.h"
#include "tool/macro.h"

#define TINKER_F1EXTN(SUFFIX, HAS_IMPL, RETURN_TYPE, FUNC, ...)                                    \
   TINKER_F1EXTN_##HAS_IMPL##__(RETURN_TYPE, FUNC, SUFFIX, __VA_ARGS__)
#define TINKER_F1EXTN_1__(R, F, S, ...) extern R F##_##S(__VA_ARGS__)
#define TINKER_F1EXTN_0__(R, F, S, ...)                                                            \
   MAYBE_UNUSED static R F##_##S(__VA_ARGS__)                                                      \
   {                                                                                               \
      throw FatalError("F1EXTN: " #R " " #F "_" #S "(...) is not implemented.");                   \
   }

#define TINKER_F1CALL(SUFFIX1, FUNC, ...) FUNC##_##SUFFIX1(__VA_ARGS__)

#define TINKER_F2EXTN(SUFFIX1, HAS_IMPL1, SUFFIX2, HAS_IMPL2, RETURN_TYPE, FUNC, ...)              \
   TINKER_F1EXTN(SUFFIX1, HAS_IMPL1, RETURN_TYPE, FUNC, __VA_ARGS__);                              \
   TINKER_F1EXTN(SUFFIX2, HAS_IMPL2, RETURN_TYPE, FUNC, __VA_ARGS__)

#define TINKER_F2PICK(SUFFIX1, SUFFIX2, FUNC, ...)                                                 \
   if (0) {                                                                                        \
   } else if (TINKER_F2PICK_##SUFFIX1##__) {                                                       \
      return FUNC##_##SUFFIX1(__VA_ARGS__);                                                        \
   } else if (TINKER_F2PICK_##SUFFIX2##__) {                                                       \
      return FUNC##_##SUFFIX2(__VA_ARGS__);                                                        \
   } else if (TINKER_HOST) {                                                                       \
      return FUNC##_acc(__VA_ARGS__);                                                              \
   } else {                                                                                        \
      throw FatalError("F2PICK: " #FUNC "(...) does not have a fallback implementation.");         \
   }
#define TINKER_F2PICK_cu__  TINKER_GPULANG_CUDA
#define TINKER_F2PICK_acc__ TINKER_GPULANG_OPENACC
