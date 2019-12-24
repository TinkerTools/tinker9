#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
namespace platform {
constexpr int cpu = 0x0001;
constexpr int gpu = 0x0100;
constexpr int acc = 0x0200;
constexpr int cuda = 0x0400;


extern int config;
}
namespace p = platform;
TINKER_NAMESPACE_END
