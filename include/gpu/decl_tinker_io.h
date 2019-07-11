#ifndef TINKER_GPU_DECL_TINKER_IO_
#define TINKER_GPU_DECL_TINKER_IO_

#include "decl_basic.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void copyin_tinker_arc(const std::string& arcfile, int first1, int last1,
                       int step);
}
TINKER_NAMESPACE_END

#endif
