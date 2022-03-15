#pragma once
#include "image3.h"

/**
 * \def image
 * \ingroup math
 * Apply periodic boundary conditions to displacement (`xr, yr, zr`) and
 * preserve the correct signs.
 *
 * For testing purpose, defining it in the source code before including this
 * header file will overwrite the default macro definition. Otherwise, it needs
 * to be undef-ed.
 */
#ifndef image
#   define image(x, y, z) image_general(x, y, z, TINKER_IMAGE_ARGS)
#endif

#ifndef image2
#   define image2(x, y, z) image2_general(x, y, z, TINKER_IMAGE_ARGS)
#endif

#ifndef imagen2
#   define imagen2(x, y, z) imagen2_general(x, y, z, TINKER_IMAGE_ARGS)
#endif
