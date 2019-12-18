#pragma once
#include "box.h"
#include "macro_void_cuda_def.h"
#include "mathfunc.h"


TINKER_NAMESPACE_BEGIN
// image


#pragma acc routine seq
__device__
inline void image_triclinic(real& restrict xr, real& restrict yr,
                            real& restrict zr, const Box* restrict pb)
{
   real fx = zr * pb->recip[0][2] + yr * pb->recip[0][1] + xr * pb->recip[0][0];
   real fy = zr * pb->recip[1][2] + yr * pb->recip[1][1];
   real fz = zr * pb->recip[2][2];
   fx = REAL_FLOOR(fx + 0.5f);
   fy = REAL_FLOOR(fy + 0.5f);
   fz = REAL_FLOOR(fz + 0.5f);
   xr -= (fz * pb->lvec[0][2] + fy * pb->lvec[0][1] + fx * pb->lvec[0][0]);
   yr -= (fz * pb->lvec[1][2] + fy * pb->lvec[1][1]);
   zr -= (fz * pb->lvec[2][2]);
}


#pragma acc routine seq
__device__
inline void image_monoclinic(real& restrict xr, real& restrict yr,
                             real& restrict zr, const Box* restrict pb)
{
   real fx = zr * pb->recip[0][2] + xr * pb->recip[0][0];
   real fy = yr * pb->recip[1][1];
   real fz = zr * pb->recip[2][2];
   fx = REAL_FLOOR(fx + 0.5f);
   fy = REAL_FLOOR(fy + 0.5f);
   fz = REAL_FLOOR(fz + 0.5f);
   xr -= (fz * pb->lvec[0][2] + fx * pb->lvec[0][0]);
   yr -= (fy * pb->lvec[1][1]);
   zr -= (fz * pb->lvec[2][2]);
}


#pragma acc routine seq
__device__
inline void image_orthogonal(real& restrict xr, real& restrict yr,
                             real& restrict zr, const Box* restrict pb)
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
inline void image_null(real& restrict, real& restrict, real& restrict,
                       const Box* restrict)
{}


#pragma acc routine seq
__device__
inline void image_general(real& restrict xr, real& restrict yr,
                          real& restrict zr, const Box* restrict pb)
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


#pragma acc routine seq
__device__
inline void imagen_triclinic(real& restrict xr, real& restrict yr,
                             real& restrict zr, const Box* restrict pb)
{
   image_triclinic(xr, yr, zr, pb);
}


#pragma acc routine seq
__device__
inline void imagen_monoclinic(real& restrict xr, real& restrict yr,
                              real& restrict zr, const Box* restrict pb)
{
   image_monoclinic(xr, yr, zr, pb);
}


#pragma acc routine seq
__device__
inline void imagen_orthogonal(real& restrict xr, real& restrict yr,
                              real& restrict zr, const Box* restrict pb)
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


#pragma acc routine seq
__device__
inline void imagen_null(real& restrict, real& restrict, real& restrict,
                        const Box* restrict)
{}


#pragma acc routine seq
__device__
inline void imagen_general(real& restrict xr, real& restrict yr,
                           real& restrict zr, const Box* restrict pb)
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
TINKER_NAMESPACE_END
