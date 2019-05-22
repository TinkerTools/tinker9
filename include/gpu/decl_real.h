#ifndef TINKER_GPU_DECL_REAL_H_
#define TINKER_GPU_DECL_REAL_H_

#include "util/macro.h"
#include <string>

/// @file
/// Type `real' was defined inside `TINKER_NAMESPACE::gpu'.

TINKER_NAMESPACE_BEGIN
namespace gpu {
#if defined(TINKER_GPU_DOUBLE)
typedef double real;
#elif defined(TINKER_GPU_SINGLE)
typedef float real;
#else
static_assert(false, "");
#endif
}
TINKER_NAMESPACE_END

#endif
