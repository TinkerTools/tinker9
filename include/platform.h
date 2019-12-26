#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
namespace platform {
constexpr int acc_pltfm = 0x001;
constexpr int cu_pltfm = 0x002;


#if TINKER_HOST
constexpr int config = acc_pltfm;
#endif


#if TINKER_CUDART
constexpr int config = acc_pltfm + cu_pltfm;
#endif
}
TINKER_NAMESPACE_END
