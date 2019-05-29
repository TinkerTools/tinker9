#ifndef TINKER_GPU_ACC_IMAGE_H_
#define TINKER_GPU_ACC_IMAGE_H_

#include "decl_box.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
#pragma acc routine seq
void image(real& __restrict__ xr, real& __restrict__ yr, real& __restrict__ zr,
           const box_st* __restrict__ pb);

#pragma acc routine seq
void imagen(real& __restrict__ xr, real& __restrict__ yr, real& __restrict__ zr,
            const box_st* __restrict__ pb);

#pragma acc routine seq
real volbox(const box_st* __restrict__ pb);
}
TINKER_NAMESPACE_END

#endif
