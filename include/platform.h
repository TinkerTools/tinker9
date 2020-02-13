#pragma once
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
namespace platform {
constexpr int NOT_SET = 0x000;
constexpr int ACC_PLTFM = 0x001;
constexpr int CU_PLTFM = 0x002;


TINKER_EXTERN int config;
}


void platform_data(rc_op);
TINKER_NAMESPACE_END
