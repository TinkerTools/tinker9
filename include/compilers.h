#pragma once
#include "macro.h"
#include <string>


TINKER_NAMESPACE_BEGIN
std::string cxx_compiler_name();
std::string acc_compiler_name();
std::string cuda_compiler_name();
TINKER_NAMESPACE_END
