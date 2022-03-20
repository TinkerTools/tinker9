#pragma once
#include "glob.accasync.h"

namespace tinker {
/// \ingroup nvidia
/// \brief Number of threads in a warp.
constexpr unsigned WARP_SIZE = 32;
/// \ingroup nvidia
/// \brief Mask for all of the lanes in a warp.
constexpr unsigned ALL_LANES = 0xFFFFFFFF;
/// \ingroup nvidia
/// \brief Default dimension of thread blocks.
// constexpr unsigned BLOCK_DIM = 64;
constexpr unsigned BLOCK_DIM = 128;
// constexpr unsigned BLOCK_DIM = 256;

constexpr int PME_BLOCKDIM = 64;
TINKER_EXTERN int ndevice;
TINKER_EXTERN int idevice;
}
