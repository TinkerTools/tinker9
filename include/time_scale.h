#pragma once
#include "macro.h"
#include <map>
#include <string>


TINKER_NAMESPACE_BEGIN
using TimeScaleConfig = std::map<std::string, int>;
const TimeScaleConfig& default_tsconfig();
TINKER_NAMESPACE_END
