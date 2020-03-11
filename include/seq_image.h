#pragma once
#include "seq_image_real3.h"


#define image_3args(x, y, z)    image_general(x, y, z, TINKER_IMAGE_ARGS)
#define image_4args(x, y, z, b) image_general(x, y, z, b)
#define image_nargs(...)                                                       \
   TINKER_GET_5TH_ARG(__VA_ARGS__, image_4args, image_3args)
/**
 * \def image
 * \ingroup macro
 * Apply periodic boundarY conditions to displacements (`xr, yr, zr`) and
 * preserve the correct signs.
 *
 * For testing purpose, define it in the source code before including this
 * header file will overwrite the default definition.
 */
#ifndef image
#   define image(...) image_nargs(__VA_ARGS__)(__VA_ARGS__)
#endif
