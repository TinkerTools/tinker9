#pragma once
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
namespace platform {
constexpr int acc_pltfm = 0x001;
constexpr int cu_pltfm = 0x002;


TINKER_EXTERN int config;
}


void platform_data(rc_op);
TINKER_NAMESPACE_END
