#pragma once
#include "tool/macro.h"

/// \page async
/// \ingroup async
///
/// This is a snippet of the OpenACC code.
/// ```cpp
/// #pragma acc parallel loop
/// for (int i = 0; i < limit; ++i) {
///    array[i] = i;
/// }
/// ```
/// By default, the CPU thread will wait/block until the parallel loop finishes.
/// If you want the CPU thread to proceed without waiting for the parallel loop
/// to finish, you may add the `async` directive,
/// ```cpp
/// #pragma acc parallel loop async(queue)
/// // or #pragma acc parallel loop async
/// for (int i = 0; i < limit; ++i) {
///    array[i] = i;
/// }
/// ```
/// where `queue` is an optional hardwired integer or an integer variable.
/// A special integer constant value `acc_async_sync` is defined in the OpenACC
/// standard that can be used in the `async` directive as a queue number,
/// to obtain an synchronous/blocking behavior. Implementations may be different
/// on different platforms though, on the CUDA platform, every OpenACC queue is
/// built on top of a CUDA stream.

namespace tinker {
/// \ingroup async
/// \brief Global handles for the GPU runtime libraries.
namespace g {
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

/// \ingroup async
TINKER_EXTERN bool use_pme_stream;
}
