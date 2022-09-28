#pragma once
#include "tool/macro.h"

/// \page async Asynchronous Queues and Streams
/// \ingroup async
///
/// This is a snippet of the OpenACC code.
/// ```cpp
/// #pragma acc parallel loop
/// for (int i = 0; i < limit; ++i) {
///    array[i] = i;
/// }
/// ```
/// By default, the host thread will wait until the parallel loop finishes.
/// If you want the host thread to proceed without waiting for the parallel loop
/// to finish, you may add \c async,
/// ```cpp
/// #pragma acc parallel loop async(queue)
/// // or #pragma acc parallel loop async
/// for (int i = 0; i < limit; ++i) {
///    array[i] = i;
/// }
/// ```
/// where \c queue is an optional hardwired integer or an integer variable.
/// A special integer constant value \c acc_async_sync is defined by the OpenACC
/// standard that can be used in the \c async directive as a queue number,
/// to achieve an synchronous/blocking behavior.
/// Implementations may be different on different platforms though,
/// on the CUDA platform every OpenACC queue is built on top of a CUDA stream.

namespace tinker {
/// Global handles for the GPU runtime libraries. \ingroup async
namespace g {
TINKER_EXTERN int q0;   ///< Default OpenACC async queue. \ingroup async
TINKER_EXTERN int q1;   ///< Default OpenACC sync queue. \ingroup async
TINKER_EXTERN int qpme; ///< OpenACC async queue for %PME. \ingroup async
}
TINKER_EXTERN bool use_pme_stream; ///< Logical flag for use of a separate CUDA stream for %PME. \ingroup async
}
