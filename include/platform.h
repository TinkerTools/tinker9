#pragma once
#include "tool/rcman.h"

namespace tinker {
/// \brief Flags for the major platforms.
/// \ingroup platform
enum class Platform
{
   UNSET = 0x000, ///< Flag for the unset platform.
   ACC = 0x001,   ///< Flag for the OpenACC platform.
   CUDA = 0x002   ///< Flag for the CUDA platform.
};
TINKER_ENABLE_ENUM_BITMASK(Platform);

void platformData(RcOp);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
TINKER_EXTERN Platform pltfm_config;
}
