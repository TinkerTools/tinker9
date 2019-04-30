#ifndef TINKER_GPU_PREC_H_
#define TINKER_GPU_PREC_H_

#include "util/macro.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
#ifdef TINKER_GPU_DOUBLE
typedef double real;
#endif

#ifdef TINKER_GPU_SINGLE
typedef float real;
#endif
}
TINKER_NAMESPACE_END

#endif
