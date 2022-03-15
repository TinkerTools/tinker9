#pragma once
#include "macro.h"

namespace tinker {
/**
 * \brief Global handles for the GPU runtime libraries.
 */
namespace g {
/**
 * \page async
 * \ingroup async
 *
 * This is a helloworld-level OpenACC code snippet.
 * ```cpp
 * #pragma acc parallel loop
 * for (int i = 0; i < limit; ++i) {
 *    array[i] = i;
 * }
 * ```
 * By default, the CPU thread will wait/block until the parallel loop finishes.
 * If you want the CPU thread to proceed without waiting for the parallel loop,
 * you may add the `async` directive,
 * ```cpp
 * #pragma acc parallel loop async(queue)
 * // or #pragma acc parallel loop async
 * for (int i = 0; i < limit; ++i) {
 *    array[i] = i;
 * }
 * ```
 * where `queue` is an optional hardwired integer or integer variable.
 * A special integer constant value `acc_async_sync` is defined in the OpenACC
 * standard that can be used in the `async` directive as a queue number,
 * to obtain an synchronous/blocking behavior. Implementations may be different
 * on different platforms though, on the CUDA platform, every OpenACC queue is
 * built on top of a CUDA stream.
 */
/// \ingroup async
/// \brief Default OpenACC async queue.
TINKER_EXTERN int q0;
/// \ingroup async
/// \brief Default OpenACC sync queue.
TINKER_EXTERN int q1;
/// \ingroup async
/// \brief OpenACC async queue for %PME.
TINKER_EXTERN int qpme;
}
TINKER_EXTERN bool use_pme_stream;
}
