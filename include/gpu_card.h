#ifndef TINKER_GPU_CARD_H_
#define TINKER_GPU_CARD_H_

#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int ndevice, idevice;

void gpu_card_data(rc_op op);
int get_grid_size(int nthreads_per_block);
int get_block_size(int shared_bytes_per_thread);
TINKER_NAMESPACE_END

#endif
