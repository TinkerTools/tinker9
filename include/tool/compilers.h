#pragma once
#include <string>

namespace tinker {
/// \addtogroup platform
/// \{
std::string cxxCompilerName();  ///< Returns the name of the C++ compiler.
std::string accCompilerName();  ///< Returns the name of the OpenACC compiler.
std::string cudaCompilerName(); ///< Returns the name of the CUDA compiler.
/// \}
}
