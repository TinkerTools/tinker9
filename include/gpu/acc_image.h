#ifndef TINKER_GPU_ACC_IMAGE_H_
#define TINKER_GPU_ACC_IMAGE_H_

#include "decl_box.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
/**
 * applys periodic boundary conditions to displacement (xr, yr, zr) and
 * preserves the correct signs
 */
#pragma acc routine seq
void image(real& __restrict__ xr, real& __restrict__ yr, real& __restrict__ zr,
           const box_t* __restrict__ pb);

/**
 * applys periodic boundary conditions to displacement (xr, yr, zr) but only
 * guarantee the lengths are correct
 */
#pragma acc routine seq
void imagen(real& __restrict__ xr, real& __restrict__ yr, real& __restrict__ zr,
            const box_t* __restrict__ pb);
}
TINKER_NAMESPACE_END

#endif
