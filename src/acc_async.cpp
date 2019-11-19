#include "async.h"
#if TINKER_CUDART && defined(_OPENACC)
#   include <openacc.h>
#endif


TINKER_NAMESPACE_BEGIN
void async_acc_queue_data(rc_op op)
{
#if TINKER_CUDART && defined(_OPENACC)
   if (op & rc_dealloc) {
      async_acc = nullptr;
   }

   if (op & rc_alloc) {
      // int handle = acc_get_default_async();
      int handle = acc_async_sync;
      async_acc = acc_get_cuda_stream(handle);
   }
#endif
}
TINKER_NAMESPACE_END
