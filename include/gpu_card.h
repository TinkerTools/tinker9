#pragma once
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
/// \ingroup nvidia
/// Number of threads in a warp.
constexpr unsigned WARP_SIZE = 32;

/// \ingroup nvidia
/// Maximum number of threads allowed in a thread block.
constexpr unsigned MAX_BLOCK_SIZE = 256;

// constexpr int BLOCK_DIM = 32;
// constexpr int BLOCK_DIM = 64;
constexpr int BLOCK_DIM = 128;

TINKER_EXTERN int ndevice, idevice;

void gpu_card_data(rc_op op);
int get_grid_size(int nthreads_per_block);
int get_block_size(int shared_bytes_per_thread);
TINKER_NAMESPACE_END
/// \defgroup nvidia Nvidia GPU
