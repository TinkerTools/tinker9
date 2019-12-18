#pragma once
#include "macro.h"
#include <string>


TINKER_NAMESPACE_BEGIN
void stopwatch_start();
void stopwatch_lap(std::string log = "");
void stopwatch_stop();
void stopwatch_reset();
TINKER_NAMESPACE_END
