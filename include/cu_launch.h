#pragma once
#ifndef __CUDACC__
#   error This header can only be included inside CUDA source files.
#endif
#include "gpu_card.h"
#include <utility>


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup nvidia
 * \brief Launch a CUDA kernel with specified thread block size.
 * \param block_size Thread block size.
 * \param nparallel  Loop limit.
 * \param k          CUDA `__globa__` kernel.
 * \param args       Argument list of the kernel.
 */
template <class K, class... Ts>
void launch_kernel2(int block_size, int nparallel, K k, Ts&&... args)
{
   int gs = (nparallel + block_size - 1) / block_size;
   k<<<gs, block_size>>>(std::forward<Ts>(args)...);
}


/**
 * \ingroup nvidia
 * \brief Launch a CUDA kernel with the a thread block size of 128.
 * \param nparallel Loop limit.
 * \param k         CUDA `__globa__` kernel.
 * \param args      Argument list of the kernel.
 */
template <class K, class... Ts>
void launch_kernel1(int nparallel, K k, Ts&&... args)
{
   const int block_size = 128;
   launch_kernel2(block_size, nparallel, k, std::forward<Ts>(args)...);
}
TINKER_NAMESPACE_END
