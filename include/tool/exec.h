#pragma once
#include <string>

namespace tinker {
/// \ingroup platform
/// Executes a command in shell and returns the output in a string.
std::string exec(const std::string& cmd);
}
