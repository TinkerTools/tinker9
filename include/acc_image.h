#ifndef TINKER_ACC_IMAGE_H_
#define TINKER_ACC_IMAGE_H_

#include "box.h"

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * apply periodic boundary conditions to displacement (@c xr, @c yr, @c zr) and
 * preserve the correct signs
 */
#pragma acc routine seq
void image(real& __restrict__ xr, real& __restrict__ yr, real& __restrict__ zr,
           const Box* __restrict__ pb);

/**
 * @brief
 * apply periodic boundary conditions to displacement (@c xr, @c yr, @c zr) but
 * only guarantee the lengths are correct
 */
#pragma acc routine seq
void imagen(real& __restrict__ xr, real& __restrict__ yr, real& __restrict__ zr,
            const Box* __restrict__ pb);
TINKER_NAMESPACE_END

#endif
