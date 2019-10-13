#pragma once
#ifndef __CUDACC__
#   error This header can only work with CUDA source files.
#endif
#include "macro.h"
#include <utility>


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup nvidia
 */
template <class K, class... Ts>
void launch_kernel2(int gs, int bs, K k, Ts&&... args)
{
   k<<<gs, bs>>>(std::forward<Ts>(args)...);
}


/**
 * \ingroup nvidia
 */
template <class K, class... Ts>
void launch_kernel3(int gs, int bs, int ss, K k, Ts&&... args)
{
   k<<<gs, bs, ss>>>(std::forward<Ts>(args)...);
}
TINKER_NAMESPACE_END
