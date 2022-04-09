#pragma once
#include "tool/error.h"
#include "tool/macro.h"

#define TINKER_F1EXTRN(RETURN_TYPE, FUNC, SUFFIX, HAS_IMPL, ...)                                   \
   TINKER_F1EXTRN_##HAS_IMPL##__(RETURN_TYPE, FUNC, SUFFIX, __VA_ARGS__)
#define TINKER_F1EXTRN_1__(R, F, S, ...) extern R F##_##S(__VA_ARGS__)
#define TINKER_F1EXTRN_0__(R, F, S, ...)                                                           \
   static R F##_##S(__VA_ARGS__)                                                                   \
   {                                                                                               \
      throw FatalError("F1EXTRN: " #F "_" #S "(...) is not implemented.");                         \
   }
#define TINKER_F2EXTRN(R, F, S1, H1, S2, H2, ...)                                                  \
   TINKER_F1EXTRN(R, F, S1, H1, __VA_ARGS__);                                                      \
   TINKER_F1EXTRN(R, F, S2, H2, __VA_ARGS__)

#define TINKER_F2CALL(F, S1, S2, ...)                                                              \
   if (0) {                                                                                        \
   } else if (TINKER_F2CALL_##S1) {                                                                \
      return F##_##S1(__VA_ARGS__);                                                                \
   } else if (TINKER_F2CALL_##S2) {                                                                \
      return F##_##S2(__VA_ARGS__);                                                                \
   } else if (TINKER_HOST) {                                                                       \
      return F##_acc(__VA_ARGS__);                                                                 \
   } else {                                                                                        \
      throw FatalError("F2CALL: " #F "(...) does not have a fallback implementation.");            \
   }
#define TINKER_F2CALL_cu  TINKER_GPULANG_CUDA
#define TINKER_F2CALL_acc TINKER_GPULANG_OPENACC
