#pragma once
#include "tool/macro.h"
#include "tool/platform.h"

namespace tinker {
void throwExceptionMissingFunction(const char* functionName, const char* file,
   int lineNum);
}

#if TINKER_GPULANG_OPENACC // mixed source code: openacc and cuda

#   define TINKER_FEXTN1_acc1__(R, F, ...) \
      TINKER_FEXTN1_NORMAL_(R, F, acc, __VA_ARGS__)
#   define TINKER_FEXTN1_acc0__(R, F, ...) TINKER_FEXTN1_EMPTY__
#   define TINKER_FEXTN1_cu1__(R, F, ...) \
      TINKER_FEXTN1_NORMAL_(R, F, cu, __VA_ARGS__)
#   define TINKER_FEXTN1_cu0__(R, F, ...) TINKER_FEXTN1_EMPTY__

#   define TINKER_FCALL0(F, ...) F##_acc(__VA_ARGS__)

#   define TINKER_FCALL2_acc0_cu1__(F, ...) \
      TINKER_FCALL0_NORMAL_(F, cu, __VA_ARGS__)
#   define TINKER_FCALL2_acc1_cu0__(F, ...) \
      TINKER_FCALL0_NORMAL_(F, acc, __VA_ARGS__)
#   define TINKER_FCALL2_acc1_cu1__(F, ...)          \
      (pltfm_config & Platform::CUDA)                \
         ? TINKER_FCALL0_NORMAL_(F, cu, __VA_ARGS__) \
         : TINKER_FCALL0_NORMAL_(F, acc, __VA_ARGS__)

#elif TINKER_GPULANG_CUDA // pure cuda

#   define TINKER_FEXTN1_acc1__(R, F, ...) TINKER_FEXTN1_EMPTY__
#   define TINKER_FEXTN1_acc0__(R, F, ...) TINKER_FEXTN1_EMPTY__
#   define TINKER_FEXTN1_cu1__(R, F, ...) \
      TINKER_FEXTN1_NORMAL_(R, F, cu, __VA_ARGS__)
#   define TINKER_FEXTN1_cu0__(R, F, ...) TINKER_FEXTN1_EMPTY__

#   define TINKER_FCALL0(F, ...) F##_cu(__VA_ARGS__)

#   define TINKER_FCALL2_acc0_cu1__(F, ...) \
      TINKER_FCALL0_NORMAL_(F, cu, __VA_ARGS__)
#   define TINKER_FCALL2_acc1_cu0__(F, ...) TINKER_FCALL0_ERROR__(F, cu)
#   define TINKER_FCALL2_acc1_cu1__(F, ...) \
      TINKER_FCALL0_NORMAL_(F, cu, __VA_ARGS__)

#else // host

#   define TINKER_FEXTN1_acc1__(R, F, ...) \
      TINKER_FEXTN1_NORMAL_(R, F, acc, __VA_ARGS__)
#   define TINKER_FEXTN1_acc0__(R, F, ...) TINKER_FEXTN1_EMPTY__
#   define TINKER_FEXTN1_cu1__(R, F, ...)  TINKER_FEXTN1_EMPTY__
#   define TINKER_FEXTN1_cu0__(R, F, ...)  TINKER_FEXTN1_EMPTY__

#   define TINKER_FCALL0(F, ...) F##_acc(__VA_ARGS__)

#   define TINKER_FCALL2_acc0_cu1__(F, ...) TINKER_FCALL0_ERROR__(F, cu)
#   define TINKER_FCALL2_acc1_cu0__(F, ...) \
      TINKER_FCALL0_NORMAL_(F, acc, __VA_ARGS__)
#   define TINKER_FCALL2_acc1_cu1__(F, ...) \
      TINKER_FCALL0_NORMAL_(F, acc, __VA_ARGS__)

#endif // end

#define TINKER_FVOID2(SUFFIX_HAS_IMPL1, SUFFIX_HAS_IMPL2, FUNC, ...) \
   TINKER_FEXTN2__(SUFFIX_HAS_IMPL1, SUFFIX_HAS_IMPL2, void, FUNC, __VA_ARGS__)

#define TINKER_FEXTN2__(SH1, SH2, R, F, ...) \
   TINKER_FEXTN1__(SH1, R, F, __VA_ARGS__);  \
   TINKER_FEXTN1__(SH2, R, F, __VA_ARGS__)

#define TINKER_FEXTN1__(SH, R, F, ...) TINKER_FEXTN1_##SH##__(R, F, __VA_ARGS__)

#define TINKER_FVOID1 TINKER_FVOID2
#define TINKER_FCALL1 TINKER_FCALL2

#define TINKER_FCALL2(SUFFIX_HAS_IMPL1, SUFFIX_HAS_IMPL2, FUNC, ...) \
   TINKER_FCALL2_##SUFFIX_HAS_IMPL1##_##SUFFIX_HAS_IMPL2##__(FUNC, __VA_ARGS__)

#define TINKER_FEXTN1_EMPTY__
#define TINKER_FEXTN1_NORMAL_(R, F, S, ...) extern R F##_##S(__VA_ARGS__)

#define TINKER_FCALL0_ERROR__(F, S) \
   throwExceptionMissingFunction(#F "_" #S, __FILE__, __LINE__)
#define TINKER_FCALL0_NORMAL_(F, S, ...) F##_##S(__VA_ARGS__)
