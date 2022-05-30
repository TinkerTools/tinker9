#pragma once
#include <string>

namespace tinker {
/// \addtogroup platform
/// \{
std::string cxxCompilerName();  ///< \return  Name of the C++ compiler.
std::string accCompilerName();  ///< \return  Name of the OpenACC compiler.
std::string cudaCompilerName(); ///< \return  Name of the CUDA compiler.
/// \}
}
