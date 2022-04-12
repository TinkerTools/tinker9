#pragma once
#include "tool/macro.h"

namespace tinker {
void throwExceptionMissingFunction(const char* functionName);
}

#if TINKER_GPULANG_OPENACC // mixed source code: openacc and cuda

#   define TINKER_FEXTN1__(S, H, R, F, ...) TINKER_FEXTN1_##H##__(R, F, S, __VA_ARGS__)
#   define TINKER_FEXTN1_1__(R, F, S, ...)  TINKER_FEXTN1_NORMAL_(R, F, S, __VA_ARGS__)
#   define TINKER_FEXTN1_0__(R, F, S, ...)  TINKER_FEXTN1_EMPTY__

#   define TINKER_FCALL0(F, ...) F##_acc(__VA_ARGS__)

#   define TINKER_FCALL2_cu0_acc1__(F, ...) TINKER_FCALL1_NORMAL_(F, acc, __VA_ARGS__)
#   define TINKER_FCALL2_cu1_acc0__(F, ...) TINKER_FCALL1_NORMAL_(F, cu, __VA_ARGS__)
#   define TINKER_FCALL2_cu1_acc1__(F, ...)                                                        \
      (pltfm_config & Platform::CUDA) ? TINKER_FCALL1_NORMAL_(F, cu, __VA_ARGS__)                  \
                                      : TINKER_FCALL1_NORMAL_(F, acc, __VA_ARGS__)

#elif TINKER_GPULANG_CUDA // pure cuda

#   define TINKER_FEXTN1__(S, H, R, F, ...) TINKER_FEXTN1_##S##H##__(R, F, __VA_ARGS__)
#   define TINKER_FEXTN1_acc1__(R, F, ...)  TINKER_FEXTN1_EMPTY__
#   define TINKER_FEXTN1_acc0__(R, F, ...)  TINKER_FEXTN1_EMPTY__
#   define TINKER_FEXTN1_cu1__(R, F, ...)   TINKER_FEXTN1_NORMAL_(R, F, cu, __VA_ARGS__)
#   define TINKER_FEXTN1_cu0__(R, F, ...)   TINKER_FEXTN1_EMPTY__

#   define TINKER_FCALL0(F, ...) F##_cu(__VA_ARGS__)

#   define TINKER_FCALL2_cu0_acc1__(F, ...) TINKER_FCALL1_ERROR__(F, cu)
#   define TINKER_FCALL2_cu1_acc0__(F, ...) TINKER_FCALL1_NORMAL_(F, cu, __VA_ARGS__)
#   define TINKER_FCALL2_cu1_acc1__(F, ...) TINKER_FCALL1_NORMAL_(F, cu, __VA_ARGS__)

#else // host

#   define TINKER_FEXTN1__(S, H, R, F, ...) TINKER_FEXTN1_##S##H##__(R, F, __VA_ARGS__)
#   define TINKER_FEXTN1_acc1__(R, F, ...)  TINKER_FEXTN1_NORMAL_(R, F, acc, __VA_ARGS__)
#   define TINKER_FEXTN1_acc0__(R, F, ...)  TINKER_FEXTN1_EMPTY__
#   define TINKER_FEXTN1_cu1__(R, F, ...)   TINKER_FEXTN1_EMPTY__
#   define TINKER_FEXTN1_cu0__(R, F, ...)   TINKER_FEXTN1_EMPTY__

#   define TINKER_FCALL0(F, ...) F##_acc(__VA_ARGS__)

#   define TINKER_FCALL2_cu1_acc0__(F, ...) TINKER_FCALL1_ERROR__(F, cu)
#   define TINKER_FCALL2_cu0_acc1__(F, ...) TINKER_FCALL1_NORMAL_(F, acc, __VA_ARGS__)
#   define TINKER_FCALL2_cu1_acc1__(F, ...) TINKER_FCALL1_NORMAL_(F, acc, __VA_ARGS__)

#endif // end

#define TINKER_FCALL2_acc0_cu1__ TINKER_FCALL2_cu1_acc0__
#define TINKER_FCALL2_acc1_cu0__ TINKER_FCALL2_cu0_acc1__
#define TINKER_FCALL2_acc1_cu1__ TINKER_FCALL2_cu1_acc1__

#define TINKER_FVOID2(SUFFIX1, HAS_IMPL1, SUFFIX2, HAS_IMPL2, FUNC, ...)                           \
   TINKER_FEXTN2__(SUFFIX1, HAS_IMPL1, SUFFIX2, HAS_IMPL2, void, FUNC, __VA_ARGS__)

#define TINKER_FCALL2(SUFFIX1, HAS_IMPL1, SUFFIX2, HAS_IMPL2, FUNC, ...)                           \
   TINKER_FCALL2_##SUFFIX1##HAS_IMPL1##_##SUFFIX2##HAS_IMPL2##__(FUNC, __VA_ARGS__)

#define TINKER_FEXTN2__(S1, H1, S2, H2, R, F, ...)                                                 \
   TINKER_FEXTN1__(S1, H1, R, F, __VA_ARGS__);                                                     \
   TINKER_FEXTN1__(S2, H2, R, F, __VA_ARGS__)

#define TINKER_FEXTN1_EMPTY__
#define TINKER_FEXTN1_NORMAL_(R, F, S, ...) extern R F##_##S(__VA_ARGS__)

#define TINKER_FCALL1_ERROR__(F, S)      throwExceptionMissingFunction(#F "_" #S)
#define TINKER_FCALL1_NORMAL_(F, S, ...) F##_##S(__VA_ARGS__)
