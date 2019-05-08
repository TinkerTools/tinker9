#ifndef TINKER_GPU_ACC_H_
#define TINKER_GPU_ACC_H_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
extern int queue_b;
extern int queue_nb;
int new_acc_async_queue();
void async_launches_begin(int* queue);
void async_launches_end(int queue);
void async_launches_end();
}
TINKER_NAMESPACE_END

#endif
