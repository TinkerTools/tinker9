#pragma once
#include "image3.h"


/**
 * \def image
 * \ingroup math
 * Apply periodic boundary conditions to displacements (`xr, yr, zr`) and
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


namespace tinker {
struct IMG_ORTHO
{
   SEQ_ROUTINE
   static real img2(real& restrict xr, real& restrict yr, real& restrict zr,
                    real3 l1, real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
   {
      return image2_orthogonal(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   }
};


struct IMG_MONO
{
   SEQ_ROUTINE
   static real img2(real& restrict xr, real& restrict yr, real& restrict zr,
                    real3 l1, real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
   {
      return image2_monoclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   }
};


struct IMG_TRI
{
   SEQ_ROUTINE
   static real img2(real& restrict xr, real& restrict yr, real& restrict zr,
                    real3 l1, real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
   {
      return image2_triclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   }
};


struct IMG_OCT
{
   SEQ_ROUTINE
   static real img2(real& restrict xr, real& restrict yr, real& restrict zr,
                    real3 l1, real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
   {
      return image2_oct(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   }
};
}
