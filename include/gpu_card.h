#pragma once
#include "rc_man.h"
#if TINKER_CUDART
#   include "gpu_card_cudart.h"
#endif


TINKER_NAMESPACE_BEGIN
/// \ingroup nvidia
/// Number of threads in a warp.
constexpr unsigned WARP_SIZE = 32;


/// \ingroup nvidia
/// Mask for all of the lanes in a warp.
constexpr unsigned ALL_LANES = 0xFFFFFFFF;


/// \ingroup nvidia
/// Default dimension of thread blocks.
// constexpr unsigned BLOCK_DIM = 64;
constexpr unsigned BLOCK_DIM = 128;
// constexpr unsigned BLOCK_DIM = 256;


TINKER_EXTERN int ndevice, idevice;


void gpu_card_data(rc_op op);
int get_grid_size(int nthreads_per_block);
int get_block_size(int shared_bytes_per_thread);
TINKER_NAMESPACE_END
