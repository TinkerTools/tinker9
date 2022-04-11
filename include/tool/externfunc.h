#pragma once
#include "tool/macro.h"

namespace tinker {
void throwExceptionMissingFunction(const char* functionName);
}

#if TINKER_GPULANG_OPENACC // mixed source code: openacc and cuda

#   define TINKER_F1EXTN__(S, H, R, F, ...) TINKER_F1EXTN_##H##__(R, F, S, __VA_ARGS__)
#   define TINKER_F1EXTN_1__(R, F, S, ...)  TINKER_F1EXTN_NORMAL_(R, F, S, __VA_ARGS__)
#   define TINKER_F1EXTN_0__(R, F, S, ...)  TINKER_F1EXTN_EMPTY__

#   define TINKER_F2CALL_cu0_acc1__(F, ...) TINKER_F1CALL_NORMAL_(F, acc, __VA_ARGS__)
#   define TINKER_F2CALL_cu1_acc0__(F, ...) TINKER_F1CALL_NORMAL_(F, cu, __VA_ARGS__)
#   define TINKER_F2CALL_cu1_acc1__(F, ...)                                                        \
      (pltfm_config & Platform::CUDA) ? TINKER_F1CALL_NORMAL_(F, cu, __VA_ARGS__)                  \
                                      : TINKER_F1CALL_NORMAL_(F, acc, __VA_ARGS__)

#elif TINKER_GPULANG_CUDA // pure cuda

#   define TINKER_F1EXTN__(S, H, R, F, ...) TINKER_F1EXTN_##S##H##__(R, F, __VA_ARGS__)
#   define TINKER_F1EXTN_acc1__(R, F, ...)  TINKER_F1EXTN_ERROR__
#   define TINKER_F1EXTN_acc0__(R, F, ...)  TINKER_F1EXTN_EMPTY__
#   define TINKER_F1EXTN_cu1__(R, F, ...)   TINKER_F1EXTN_NORMAL_(R, F, cu, __VA_ARGS__)
#   define TINKER_F1EXTN_cu0__(R, F, ...)   TINKER_F1EXTN_ERROR__

#   define TINKER_F2CALL_cu0_acc1__(F, ...) TINKER_F1CALL_ERROR__(F, cu)
#   define TINKER_F2CALL_cu1_acc0__(F, ...) TINKER_F1CALL_NORMAL_(F, cu, __VA_ARGS__)
#   define TINKER_F2CALL_cu1_acc1__(F, ...) TINKER_F1CALL_NORMAL_(F, cu, __VA_ARGS__)

#else // host

#   define TINKER_F1EXTN__(S, H, R, F, ...) TINKER_F1EXTN_##S##H##__(R, F, __VA_ARGS__)
#   define TINKER_F1EXTN_acc1__(R, F, ...)  TINKER_F1EXTN_NORMAL_(R, F, acc, __VA_ARGS__)
#   define TINKER_F1EXTN_acc0__(R, F, ...)  TINKER_F1EXTN_EMPTY__
#   define TINKER_F1EXTN_cu1__(R, F, ...)   TINKER_F1EXTN_ERROR__
#   define TINKER_F1EXTN_cu0__(R, F, ...)   TINKER_F1EXTN_EMPTY__

#   define TINKER_F2CALL_cu1_acc0__(F, ...) TINKER_F1CALL_ERROR__(F, cu)
#   define TINKER_F2CALL_cu0_acc1__(F, ...) TINKER_F1CALL_NORMAL_(F, acc, __VA_ARGS__)
#   define TINKER_F2CALL_cu1_acc1__(F, ...) TINKER_F1CALL_NORMAL_(F, acc, __VA_ARGS__)

#endif // end

#define TINKER_F2CALL_acc0_cu1__ TINKER_F2CALL_cu1_acc0__
#define TINKER_F2CALL_acc1_cu0__ TINKER_F2CALL_cu0_acc1__
#define TINKER_F2CALL_acc1_cu1__ TINKER_F2CALL_cu1_acc1__

#define TINKER_F2VOID(SUFFIX1, HAS_IMPL1, SUFFIX2, HAS_IMPL2, FUNC, ...)                           \
   TINKER_F2EXTN__(SUFFIX1, HAS_IMPL1, SUFFIX2, HAS_IMPL2, void, FUNC, __VA_ARGS__)

#define TINKER_F2CALL(SUFFIX1, HAS_IMPL1, SUFFIX2, HAS_IMPL2, FUNC, ...)                           \
   TINKER_F2CALL_##SUFFIX1##HAS_IMPL1##_##SUFFIX2##HAS_IMPL2##__(FUNC, __VA_ARGS__)

#define TINKER_F2EXTN__(S1, H1, S2, H2, R, F, ...)                                                 \
   TINKER_F1EXTN__(S1, H1, R, F, __VA_ARGS__);                                                     \
   TINKER_F1EXTN__(S2, H2, R, F, __VA_ARGS__)

#define TINKER_F1EXTN_EMPTY__
#define TINKER_F1EXTN_NORMAL_(R, F, S, ...) extern R F##_##S(__VA_ARGS__)
#define TINKER_F1EXTN_ERROR__               extern void throwExceptionMissingFunction(const char*)

#define TINKER_F1CALL_ERROR__(F, S)      throwExceptionMissingFunction(#F "_" #S)
#define TINKER_F1CALL_NORMAL_(F, S, ...) F##_##S(__VA_ARGS__)
