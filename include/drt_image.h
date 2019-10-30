#pragma once
#include "box.h"
#include "macro_void_cuda_def.h"
#include "mathfunc.h"


TINKER_NAMESPACE_BEGIN
#pragma acc routine seq
__device__
inline void image_null(real& RESTRICT, real& RESTRICT, real& RESTRICT,
                       const Box* RESTRICT)
{}


#pragma acc routine seq
__device__
inline void imagen_null(real& RESTRICT, real& RESTRICT, real& RESTRICT,
                        const Box* RESTRICT)
{}


// orthogonal


#pragma acc routine seq
__device__
inline void frac_orthogonal(real& RESTRICT fx, real& RESTRICT fy,
                            real& RESTRICT fz, real xr, real yr, real zr,
                            const Box* RESTRICT pb)
{
   fx = xr * pb->recip[0][0];
   fy = yr * pb->recip[1][1];
   fz = zr * pb->recip[2][2];
   fx -= REAL_FLOOR(fx + 0.5f);
   fy -= REAL_FLOOR(fy + 0.5f);
   fz -= REAL_FLOOR(fz + 0.5f);
}


#pragma acc routine seq
__device__
inline void frac_image_orthogonal(real& restrict fx, real& restrict fy,
                                  real& restrict fz, const Box* restrict pb)
{
   fx -= REAL_FLOOR(fx + 0.5f);
   fy -= REAL_FLOOR(fy + 0.5f);
   fz -= REAL_FLOOR(fz + 0.5f);
   fx *= pb->lvec[0][0];
   fy *= pb->lvec[1][1];
   fz *= pb->lvec[2][2];
}


#pragma acc routine seq
__device__
inline void image_orthogonal(real& RESTRICT xr, real& RESTRICT yr,
                             real& RESTRICT zr, const Box* RESTRICT pb)
{
   real fx = xr * pb->recip[0][0];
   real fy = yr * pb->recip[1][1];
   real fz = zr * pb->recip[2][2];
   fx = REAL_FLOOR(fx + 0.5f);
   fy = REAL_FLOOR(fy + 0.5f);
   fz = REAL_FLOOR(fz + 0.5f);
   xr -= (fx * pb->lvec[0][0]);
   yr -= (fy * pb->lvec[1][1]);
   zr -= (fz * pb->lvec[2][2]);
}


#pragma acc routine seq
__device__
inline void imagen_orthogonal(real& RESTRICT xr, real& RESTRICT yr,
                              real& RESTRICT zr, const Box* RESTRICT pb)
{
   real lx = pb->lvec[0][0];
   real lx2 = lx * 0.5f;
   real a = REAL_ABS(xr);
   xr = (a > lx2 ? a - lx : a);

   real ly = pb->lvec[1][1];
   real ly2 = ly * 0.5f;
   real b = REAL_ABS(yr);
   yr = (b > ly2 ? b - ly : b);

   real lz = pb->lvec[2][2];
   real lz2 = lz * 0.5f;
   real c = REAL_ABS(zr);
   zr = (c > lz2 ? c - lz : c);
}


// monoclinic


#pragma acc routine seq
__device__
inline void frac_monoclinic(real& RESTRICT fx, real& RESTRICT fy,
                            real& RESTRICT fz, real xr, real yr, real zr,
                            const Box* RESTRICT pb)
{
   fx = xr * pb->recip[0][0] + zr * pb->recip[0][2];
   fy = yr * pb->recip[1][1];
   fz = zr * pb->recip[2][2];
   fx -= REAL_FLOOR(fx + 0.5f);
   fy -= REAL_FLOOR(fy + 0.5f);
   fz -= REAL_FLOOR(fz + 0.5f);
}


#pragma acc routine seq
__device__
inline void frac_image_monoclinic(real& restrict fx, real& restrict fy,
                                  real& restrict fz, const Box* restrict pb)
{
   fx -= REAL_FLOOR(fx + 0.5f);
   fy -= REAL_FLOOR(fy + 0.5f);
   fz -= REAL_FLOOR(fz + 0.5f);
   fx = fx * pb->lvec[0][0] + fz * pb->lvec[0][2];
   fy *= pb->lvec[1][1];
   fz *= pb->lvec[2][2];
}


#pragma acc routine seq
__device__
inline void image_monoclinic(real& RESTRICT xr, real& RESTRICT yr,
                             real& RESTRICT zr, const Box* RESTRICT pb)
{
   real fx = xr * pb->recip[0][0] + zr * pb->recip[0][2];
   real fy = yr * pb->recip[1][1];
   real fz = zr * pb->recip[2][2];
   fx = REAL_FLOOR(fx + 0.5f);
   fy = REAL_FLOOR(fy + 0.5f);
   fz = REAL_FLOOR(fz + 0.5f);
   xr -= (fx * pb->lvec[0][0] + fz * pb->lvec[0][2]);
   yr -= (fy * pb->lvec[1][1]);
   zr -= (fz * pb->lvec[2][2]);
}


#pragma acc routine seq
__device__
inline void imagen_monoclinic(real& RESTRICT xr, real& RESTRICT yr,
                              real& RESTRICT zr, const Box* RESTRICT pb)
{
   // TODO: a real imagen routine
   image_monoclinic(xr, yr, zr, pb);
}


// triclinic


#pragma acc routine seq
__device__
inline void frac_triclinic(real& RESTRICT fx, real& RESTRICT fy,
                           real& RESTRICT fz, real xr, real yr, real zr,
                           const Box* RESTRICT pb)
{
   fx = xr * pb->recip[0][0] + yr * pb->recip[0][1] + zr * pb->recip[0][2];
   fy = yr * pb->recip[1][1] + zr * pb->recip[1][2];
   fz = zr * pb->recip[2][2];
   fx -= REAL_FLOOR(fx + 0.5f);
   fy -= REAL_FLOOR(fy + 0.5f);
   fz -= REAL_FLOOR(fz + 0.5f);
}


#pragma acc routine seq
__device__
inline void frac_image_triclinic(real& restrict fx, real& restrict fy,
                                 real& restrict fz, const Box* restrict pb)
{
   fx -= REAL_FLOOR(fx + 0.5f);
   fy -= REAL_FLOOR(fy + 0.5f);
   fz -= REAL_FLOOR(fz + 0.5f);
   fx = fx * pb->lvec[0][0] + fy * pb->lvec[0][1] + fz * pb->lvec[0][2];
   fy = fy * pb->lvec[1][1] + fz * pb->lvec[1][2];
   fz *= pb->lvec[2][2];
}


#pragma acc routine seq
__device__
inline void image_triclinic(real& RESTRICT xr, real& RESTRICT yr,
                            real& RESTRICT zr, const Box* RESTRICT pb)
{
   real fx = xr * pb->recip[0][0] + yr * pb->recip[0][1] + zr * pb->recip[0][2];
   real fy = yr * pb->recip[1][1] + zr * pb->recip[1][2];
   real fz = zr * pb->recip[2][2];
   fx = REAL_FLOOR(fx + 0.5f);
   fy = REAL_FLOOR(fy + 0.5f);
   fz = REAL_FLOOR(fz + 0.5f);
   xr -= (fx * pb->lvec[0][0] + fy * pb->lvec[0][1] + fz * pb->lvec[0][2]);
   yr -= (fy * pb->lvec[1][1] + fz * pb->lvec[1][2]);
   zr -= (fz * pb->lvec[2][2]);
}


#pragma acc routine seq
__device__
inline void imagen_triclinic(real& RESTRICT xr, real& RESTRICT yr,
                             real& RESTRICT zr, const Box* RESTRICT pb)
{
   // TODO: a real imagen routine
   image_triclinic(xr, yr, zr, pb);
}


// general


#pragma acc routine seq
__device__
inline void frac_general(real& RESTRICT fx, real& RESTRICT fy,
                         real& RESTRICT fz, real xr, real yr, real zr,
                         const Box* RESTRICT pb)
{
   switch (pb->shape) {
   case Box::ortho:
      frac_orthogonal(fx, fy, fz, xr, yr, zr, pb);
      break;
   case Box::mono:
      frac_monoclinic(fx, fy, fz, xr, yr, zr, pb);
      break;
   default:
      // Box::tri
      frac_triclinic(fx, fy, fz, xr, yr, zr, pb);
      break;
   }
}


#pragma acc routine seq
__device__
inline void frac_image_general(real& restrict fx, real& restrict fy,
                               real& restrict fz, const Box* restrict pb)
{
   switch (pb->shape) {
   case Box::ortho:
      frac_image_orthogonal(fx, fy, fz, pb);
      break;
   case Box::mono:
      frac_image_monoclinic(fx, fy, fz, pb);
      break;
   default:
      // Box::tri
      frac_image_triclinic(fx, fy, fz, pb);
      break;
   }
}


#pragma acc routine seq
__device__
inline void image_general(real& RESTRICT xr, real& RESTRICT yr,
                          real& RESTRICT zr, const Box* RESTRICT pb)
{
   switch (pb->shape) {
   case Box::ortho:
      image_orthogonal(xr, yr, zr, pb);
      break;
   case Box::mono:
      image_monoclinic(xr, yr, zr, pb);
      break;
   case Box::tri:
      image_triclinic(xr, yr, zr, pb);
      break;
   default:
      // Box::null
      image_null(xr, yr, zr, pb);
      break;
   }
}


#pragma acc routine seq
__device__
inline void imagen_general(real& RESTRICT xr, real& RESTRICT yr,
                           real& RESTRICT zr, const Box* RESTRICT pb)
{
   switch (pb->shape) {
   case Box::ortho:
      imagen_orthogonal(xr, yr, zr, pb);
      break;
   case Box::mono:
      imagen_monoclinic(xr, yr, zr, pb);
      break;
   case Box::tri:
      imagen_triclinic(xr, yr, zr, pb);
      break;
   default:
      // Box::null
      imagen_null(xr, yr, zr, pb);
      break;
   }
}


/**
 * \def frac
 * \ingroup macro
 * Calculate the fractional coordinates of (`xr, yr, zr`). The range of the
 * fractional coordinate is `[-1/2, 1/2)`.
 */
#ifndef frac
#   define frac frac_general
#endif


/**
 * \def frac_image
 * \ingroup macro
 */
#ifndef frac_image
#   define frac_image frac_image_general
#endif


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
