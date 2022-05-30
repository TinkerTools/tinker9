#pragma once
#include "tool/rcman.h"

namespace tinker {
/// \ingroup platform
/// Flags for the major platforms.
enum class Platform
{
   UNKNOWN = 0x000, ///< Flag for the unknown platform.
   ACC = 0x001,     ///< Flag for the OpenACC platform.
   CUDA = 0x002     ///< Flag for the CUDA platform.
};
TINKER_ENABLE_ENUM_BITMASK(Platform);

/// \ingroup platform
void platformData(RcOp);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
/// \ingroup platform
TINKER_EXTERN Platform pltfm_config; ///< Platform in use.
}
