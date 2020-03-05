#pragma once
#include "box.h"
#include "mathfunc.h"
#include "seq_def.h"


TINKER_NAMESPACE_BEGIN
// image


SEQ_ROUTINE
inline void image_triclinic(real& restrict xr, real& restrict yr,
                            real& restrict zr, real3 l1, real3 l2, real3 l3,
                            real3 ra, real3 rb, real3 rc)
{
   real fx = REAL_FLOOR(0.5f + zr * ra.z + yr * ra.y + xr * ra.x);
   real fy = REAL_FLOOR(0.5f + zr * rb.z + yr * rb.y);
   real fz = REAL_FLOOR(0.5f + zr * rc.z);
   xr -= (fz * l1.z + fy * l1.y + fx * l1.x);
   yr -= (fz * l2.z + fy * l2.y);
   zr -= (fz * l3.z);
}


SEQ_ROUTINE
inline void image_monoclinic(real& restrict xr, real& restrict yr,
                             real& restrict zr, real3 l1, real3 l2, real3 l3,
                             real3 ra, real3 rb, real3 rc)
{
   real fx = REAL_FLOOR(0.5f + zr * ra.z + xr * ra.x);
   real fy = REAL_FLOOR(0.5f + yr * rb.y);
   real fz = REAL_FLOOR(0.5f + zr * rc.z);
   xr -= (fz * l1.z + fx * l1.x);
   yr -= (fy * l2.y);
   zr -= (fz * l3.z);
}


SEQ_ROUTINE
inline void image_orthogonal(real& restrict xr, real& restrict yr,
                             real& restrict zr, real3 l1, real3 l2, real3 l3,
                             real3 ra, real3 rb, real3 rc)
{
   real fx = REAL_FLOOR(0.5f + xr * ra.x);
   real fy = REAL_FLOOR(0.5f + yr * rb.y);
   real fz = REAL_FLOOR(0.5f + zr * rc.z);
   xr -= (fx * l1.x);
   yr -= (fy * l2.y);
   zr -= (fz * l3.z);
}


SEQ_ROUTINE
inline void image_general(real& restrict xr, real& restrict yr,
                          real& restrict zr, real3 l1, real3 l2, real3 l3,
                          real3 ra, real3 rb, real3 rc)
{
   if (l1.z == 0) {
      image_orthogonal(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (l1.y == 0) {
      image_monoclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else {
      image_triclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   }
}


#if 0
/**
 * This imagen routine is very cunning: it is faster than the orthogonal image
 * routine but only works if two atom images are not terribly far so that their
 * distance can be corrected by 1x box size. Therefore, if we have to use this
 * routine, we either have to:
 *    - call bounds() (which has not been implemented on GPU) for every step;
 *    - or call bounds() on GPU on a regular basis and keep our fingers crossed
 *    that atoms don't diffuse too far from the PBC box.
 */
SEQ_ROUTINE
inline void imagen_orthogonal(real& xr, real& yr, real& zr, real3 l1, real3 l2,
                              real3 l3)
{
   real lx = l1.x;
   real lx2 = lx * 0.5f;
   real a = REAL_ABS(xr);
   xr = (a > lx2 ? a - lx : a);

   real ly = l2.y;
   real ly2 = ly * 0.5f;
   real b = REAL_ABS(yr);
   yr = (b > ly2 ? b - ly : b);

   real lz = l3.z;
   real lz2 = lz * 0.5f;
   real c = REAL_ABS(zr);
   zr = (c > lz2 ? c - lz : c);
}
#endif
TINKER_NAMESPACE_END


TINKER_NAMESPACE_BEGIN
// image2
#ifndef image2
#   define image2(x, y, z) image2_general(x, y, z, TINKER_IMAGE_ARGS)
#endif


SEQ_ROUTINE
inline real image2_general(real& restrict xr, real& restrict yr,
                           real& restrict zr, real3 l1, real3 l2, real3 l3,
                           real3 ra, real3 rb, real3 rc)
{
   if (l1.z == 0) {
      image_orthogonal(xr, yr, zr, l1, l2, l3, ra, rb, rc);
      return xr * xr + yr * yr + zr * zr;
   } else if (l1.y == 0) {
      image_monoclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
      return xr * xr + yr * yr + zr * zr;
   } else {
      image_triclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
      return xr * xr + yr * yr + zr * zr;
   }
}


// imagen2
#ifndef imagen2
#   define imagen2(x, y, z) imagen2_general(x, y, z, TINKER_IMAGE_ARGS)
#endif


SEQ_ROUTINE
inline real imagen2_general(real& xr, real& yr, real& zr, real3 l1, real3 l2,
                            real3 l3, real3 ra, real3 rb, real3 rc)
{
   if (l1.z == 0) {
      image_orthogonal(xr, yr, zr, l1, l2, l3, ra, rb, rc);
      return xr * xr + yr * yr + zr * zr;
   } else if (l1.y == 0) {
      image_monoclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
      return xr * xr + yr * yr + zr * zr;
   } else {
      image_triclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
      return xr * xr + yr * yr + zr * zr;
   }
}


// imagec
#ifndef imagec
#   define imagec(x0, y0, z0, xc, yc, zc)                                      \
      imagec_general(x0, y0, z0, xc, yc, zc, TINKER_IMAGE_ARGS)
#endif


SEQ_ROUTINE
inline void imagec_triclinic(real& restrict x0, real& restrict y0,
                             real& restrict z0, real xc, real yc, real zc,
                             real3 l1, real3 l2, real3 l3, real3 ra, real3 rb,
                             real3 rc)
{
   real xr = x0 - xc;
   real yr = y0 - yc;
   real zr = z0 - zc;
   real fx = REAL_FLOOR(0.5f + zr * ra.z + yr * ra.y + xr * ra.x);
   real fy = REAL_FLOOR(0.5f + zr * rb.z + yr * rb.y);
   real fz = REAL_FLOOR(0.5f + zr * rc.z);
   x0 -= (fz * l1.z + fy * l1.y + fx * l1.x);
   y0 -= (fz * l2.z + fy * l2.y);
   z0 -= (fz * l3.z);
}


SEQ_ROUTINE
inline void imagec_monoclinic(real& restrict x0, real& restrict y0,
                              real& restrict z0, real xc, real yc, real zc,
                              real3 l1, real3 l2, real3 l3, real3 ra, real3 rb,
                              real3 rc)
{
   real xr = x0 - xc;
   real yr = y0 - yc;
   real zr = z0 - zc;
   real fx = REAL_FLOOR(0.5f + zr * ra.z + xr * ra.x);
   real fy = REAL_FLOOR(0.5f + yr * rb.y);
   real fz = REAL_FLOOR(0.5f + zr * rc.z);
   x0 -= (fz * l1.z + fx * l1.x);
   y0 -= (fy * l2.y);
   z0 -= (fz * l3.z);
}


SEQ_ROUTINE
inline void imagec_orthogonal(real& restrict x0, real& restrict y0,
                              real& restrict z0, real xc, real yc, real zc,
                              real3 l1, real3 l2, real3 l3, real3 ra, real3 rb,
                              real3 rc)
{
   x0 -= REAL_FLOOR(0.5f + (x0 - xc) * ra.x) * l1.x;
   y0 -= REAL_FLOOR(0.5f + (y0 - yc) * rb.y) * l2.y;
   z0 -= REAL_FLOOR(0.5f + (z0 - zc) * rc.z) * l3.z;
}


SEQ_ROUTINE
inline void imagec_general(real& restrict x0, real& restrict y0,
                           real& restrict z0, real xc, real yc, real zc,
                           real3 l1, real3 l2, real3 l3, real3 ra, real3 rb,
                           real3 rc)
{
   if (l1.z == 0) {
      imagec_orthogonal(x0, y0, z0, xc, yc, zc, l1, l2, l3, ra, rb, rc);
   } else if (l1.y == 0) {
      imagec_monoclinic(x0, y0, z0, xc, yc, zc, l1, l2, l3, ra, rb, rc);
   } else {
      imagec_triclinic(x0, y0, z0, xc, yc, zc, l1, l2, l3, ra, rb, rc);
   }
}
TINKER_NAMESPACE_END
