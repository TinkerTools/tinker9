#ifndef TINKER_ACC_IMAGE_H_
#define TINKER_ACC_IMAGE_H_

#include "acc_mathfunc.h"
#include "box.h"

TINKER_NAMESPACE_BEGIN
#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void image_null (real& RESTRICT, real& RESTRICT, real& RESTRICT,
                        const Box* RESTRICT)
{}

#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void imagen_null (real& RESTRICT, real& RESTRICT, real& RESTRICT,
                         const Box* RESTRICT)
{}

#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void image_orthogonal (real& RESTRICT xr, real& RESTRICT yr,
                              real& RESTRICT zr, const Box* RESTRICT pb)
{
   real fx = xr * pb->recip[0][0];
   real fy = yr * pb->recip[1][1];
   real fz = zr * pb->recip[2][2];
   fx = REAL_FLOOR (fx + 0.5f);
   fy = REAL_FLOOR (fy + 0.5f);
   fz = REAL_FLOOR (fz + 0.5f);
   xr -= (fx * pb->lvec[0][0]);
   yr -= (fy * pb->lvec[1][1]);
   zr -= (fz * pb->lvec[2][2]);
}

#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void imagen_orthogonal (real& RESTRICT xr, real& RESTRICT yr,
                               real& RESTRICT zr, const Box* RESTRICT pb)
{
   real lx = pb->lvec[0][0];
   real lx2 = lx * 0.5f;
   real a = REAL_ABS (xr);
   xr = (a > lx2 ? a - lx : a);

   real ly = pb->lvec[1][1];
   real ly2 = ly * 0.5f;
   real b = REAL_ABS (yr);
   yr = (b > ly2 ? b - ly : b);

   real lz = pb->lvec[2][2];
   real lz2 = lz * 0.5f;
   real c = REAL_ABS (zr);
   zr = (c > lz2 ? c - lz : c);
}

#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void image_monoclinic (real& RESTRICT xr, real& RESTRICT yr,
                              real& RESTRICT zr, const Box* RESTRICT pb)
{
   real fx = xr * pb->recip[0][0] + zr * pb->recip[0][2];
   real fy = yr * pb->recip[1][1];
   real fz = zr * pb->recip[2][2];
   fx = REAL_FLOOR (fx + 0.5f);
   fy = REAL_FLOOR (fy + 0.5f);
   fz = REAL_FLOOR (fz + 0.5f);
   xr -= (fx * pb->lvec[0][0] + fz * pb->lvec[0][2]);
   yr -= (fy * pb->lvec[1][1]);
   zr -= (fz * pb->lvec[2][2]);
}

#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void imagen_monoclinic (real& RESTRICT xr, real& RESTRICT yr,
                               real& RESTRICT zr, const Box* RESTRICT pb)
{
   // TODO: a real imagen
   // routine
   image_monoclinic (xr, yr, zr, pb);
}

#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void image_triclinic (real& RESTRICT xr, real& RESTRICT yr,
                             real& RESTRICT zr, const Box* RESTRICT pb)
{
   real fx = xr * pb->recip[0][0] + yr * pb->recip[0][1] + zr * pb->recip[0][2];
   real fy = yr * pb->recip[1][1] + zr * pb->recip[1][2];
   real fz = zr * pb->recip[2][2];
   fx = REAL_FLOOR (fx + 0.5f);
   fy = REAL_FLOOR (fy + 0.5f);
   fz = REAL_FLOOR (fz + 0.5f);
   xr -= (fx * pb->lvec[0][0] + fy * pb->lvec[0][1] + fz * pb->lvec[0][2]);
   yr -= (fy * pb->lvec[1][1] + fz * pb->lvec[1][2]);
   zr -= (fz * pb->lvec[2][2]);
}

#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void imagen_triclinic (real& RESTRICT xr, real& RESTRICT yr,
                              real& RESTRICT zr, const Box* RESTRICT pb)
{
   // TODO: a real imagen
   // routine
   image_triclinic (xr, yr, zr, pb);
}

#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void image_general (real& RESTRICT xr, real& RESTRICT yr,
                           real& RESTRICT zr, const Box* RESTRICT pb)
{
   switch (pb->shape) {
   case Box::ortho:
      image_orthogonal (xr, yr, zr, pb);
      break;
   case Box::mono:
      image_monoclinic (xr, yr, zr, pb);
      break;
   case Box::tri:
      image_triclinic (xr, yr, zr, pb);
      break;
   default:
      // Box::null
      image_null (xr, yr, zr, pb);
      break;
   }
}

#pragma acc routine seq
CUDA_DEVICE_FUNCTION
inline void imagen_general (real& RESTRICT xr, real& RESTRICT yr,
                            real& RESTRICT zr, const Box* RESTRICT pb)
{
   switch (pb->shape) {
   case Box::ortho:
      imagen_orthogonal (xr, yr, zr, pb);
      break;
   case Box::mono:
      imagen_monoclinic (xr, yr, zr, pb);
      break;
   case Box::tri:
      imagen_triclinic (xr, yr, zr, pb);
      break;
   default:
      // Box::null
      imagen_null (xr, yr, zr, pb);
      break;
   }
}

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
#   define image image_general
#endif

/**
 * \def imagen
 * \ingroup macro
 * Apply periodic boundary conditions to displacements (`xr, yr, zr`) but only
 * guarantee the lengths are correct.
 *
 * For testing purpose, define it in the source code before including this
 * header file will overwrite the default definition.
 */
#ifndef imagen
#   define imagen imagen_general
#endif
TINKER_NAMESPACE_END

#endif
