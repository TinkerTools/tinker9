#pragma once
#include "mod/energi.h"
#include "mod/gpucard.h"
#include "tool/cudalib.h"
#include "tool/error.h"

namespace tinker {
/**
 * \ingroup cuda_syntax
 * Launch a non-blocking CUDA kernel.
 * \param st  CUDA stream.
 * \param sh  Dynamic shared memory (bytes per block).
 * \param bs  Thread block size.
 * \param np  Hint to determine the grid dimension (often is the loop limit).
 * \param k   CUDA `__global__` kernel.
 * \param a   Arguments of the kernel.
 * \see launch_k4b
 */
template <class K, class... Ts>
void launch_k4(cudaStream_t st, size_t sh, int bs, int np, K k, Ts&&... a)
{
   int gs = (np + bs - 1) / bs;
   k<<<gs, bs, sh, st>>>(std::forward<Ts>(a)...);
}

/**
 * \ingroup cuda_syntax
 * Launch a non-blocking CUDA kernel with energy, virial, and/or count buffer
 * parameters.
 * \param st  CUDA stream.
 * \param sh  Dynamic shared memory (bytes per block).
 * \param bs  Thread block size.
 * \param np  Hint to determine the grid dimension (often is the loop limit).
 * \param k   CUDA `__global__` kernel.
 * \param a   Arguments of the kernel.
 * \see launch_k4
 */
template <class K, class... Ts>
void launch_k4b(cudaStream_t st, size_t sh, int bs, int np, K k, Ts&&... a)
{
   np = std::min(np, nelem_buffer);
   int gs = (np + bs - 1) / bs;
   k<<<gs, bs, sh, st>>>(std::forward<Ts>(a)...);
}

/**
 * \ingroup cuda_syntax
 * Launch a non-blocking CUDA kernel with 0 dynamic shared memory.
 * \see launch_k4
 */
template <class K, class... Ts>
void launch_k2s(cudaStream_t st, int bs, int np, K k, Ts&&... a)
{
   const size_t sh = 0;
   launch_k4(st, sh, bs, np, k, std::forward<Ts>(a)...);
}

/**
 * \ingroup cuda_syntax
 * Launch a non-blocking CUDA kernel with 0 dynamic shared memory.
 * \see launch_k4b
 */
template <class K, class... Ts>
void launch_k2b(cudaStream_t st, int bs, int np, K k, Ts&&... a)
{
   const size_t sh = 0;
   launch_k4b(st, sh, bs, np, k, std::forward<Ts>(a)...);
}

/**
 * \ingroup cuda_syntax
 * Launch a non-blocking CUDA kernel with 0 dynamic shared memory and the
 * default block size (#BLOCK_DIM).
 * \see BLOCK_DIM, launch_k2s, launch_k4
 */
template <class K, class... Ts>
void launch_k1s(cudaStream_t st, int np, K k, Ts&&... a)
{
   const int bs = BLOCK_DIM;
   launch_k2s(st, bs, np, k, std::forward<Ts>(a)...);
}

/**
 * \ingroup cuda_syntax
 * Launch a non-blocking CUDA kernel with 0 dynamic shared memory and the
 * default block size (#BLOCK_DIM).
 * \see BLOCK_DIM, launch_k2b, launch_k4b
 */
template <class K, class... Ts>
void launch_k1b(cudaStream_t st, int np, K k, Ts&&... a)
{
   const int bs = BLOCK_DIM;
   launch_k2b(st, bs, np, k, std::forward<Ts>(a)...);
}
}
