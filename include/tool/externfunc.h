#pragma once
#include "tool/error.h"

#define TINKER_F1EXTRN(RETURN_TYPE, FUNC, SUFFIX, HAS_IMPL, ...)                                   \
   TINKER_F1EXTRN##_##HAS_IMPL##__(RETURN_TYPE, FUNC, SUFFIX, __VA_ARGS__)
#define TINKER_F1EXTRN_1__(R, F, S, ...) extern R F##_##S(__VA_ARGS__)
#define TINKER_F1EXTRN_0__(R, F, S, ...)                                                           \
   static R F##_##S(__VA_ARGS__)                                                                   \
   {                                                                                               \
      throw FatalError(#F "_" #S " is not implemented.");                                          \
      return static_cast<R>(0);                                                                    \
   }
#define TINKER_F2EXTRN(R, F, S1, H1, S2, H2, ...)                                                  \
   TINKER_F1EXTRN(R, F, S1, H1, __VA_ARGS__);                                                      \
   TINKER_F1EXTRN(R, F, S2, H2, __VA_ARGS__)
