#ifndef TINKER_GPU_DECL_BOX_H_
#define TINKER_GPU_DECL_BOX_H_

#include "mod_box.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void box_data(rc_t rc);
void box_data_copyout(const box_t& b);
}
TINKER_NAMESPACE_END

#endif
