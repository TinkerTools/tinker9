#include "gpu/acc.h"
#include <openacc.h>

TINKER_NAMESPACE_BEGIN
namespace gpu {
int queue_b;
int queue_nb;

int new_acc_async_queue() {
  static int default_ = acc_get_default_async();
  return ++default_;
}

void async_launches_begin(int* queue) {
  m_tinker_using_namespace;
  *queue = gpu::new_acc_async_queue();
}

void async_launches_end(int queue) {
  #pragma acc wait(queue)
}

void async_launches_end() {
  #pragma acc wait
}
}
TINKER_NAMESPACE_END
