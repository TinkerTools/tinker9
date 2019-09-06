#ifndef TINKER_CU_LAUNCH_H_
#define TINKER_CU_LAUNCH_H_

#include "gpu_card.h"

TINKER_NAMESPACE_BEGIN
template <class F, class... Args>
void launch_parallel_kernel(int shared_bytes_per_thread, F&& kernel,
                            Args&&... args) {
  int nthreads_per_block = get_block_size(shared_bytes_per_thread);
  int nblocks_per_grid = get_grid_size(nthreads_per_block);
  int nB_per_block = shared_bytes_per_thread * nthreads_per_block;
  kernel<<<nblocks_per_grid, nthreads_per_block, nB_per_block>>>(args...);
}
TINKER_NAMESPACE_END

#endif
