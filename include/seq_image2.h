#pragma once
#include "box.h"
#include "mathfunc.h"


TINKER_NAMESPACE_BEGIN
#pragma acc routine seq
__device__
inline real4 image2_triclinic(real xr, real yr, real zr, real3 l1, real3 l2,
                              real3 l3, real3 ra, real3 rb, real3 rc)
{
   real fx = REAL_FLOOR(0.5f + xr * ra.x + yr * ra.y + zr * ra.z);
   real fy = REAL_FLOOR(0.5f + yr * rb.y + zr * rb.z);
   real fz = REAL_FLOOR(0.5f + zr * rc.z);
   xr -= (fx * l1.x + fy * l1.y + fz * l1.z);
   yr -= (fy * l2.y + fz * l2.z);
   zr -= (fz * l3.z);
   return make_real4(xr, yr, zr, xr * xr + yr * yr + zr * zr);
}


#pragma acc routine seq
__device__
inline real4 image2_monoclinic(real xr, real yr, real zr, real3 l1, real3 l2,
                               real3 l3, real3 ra, real3 rb, real3 rc)
{
   real fx = REAL_FLOOR(0.5f + xr * ra.x + zr * ra.z);
   real fy = REAL_FLOOR(0.5f + yr * rb.y);
   real fz = REAL_FLOOR(0.5f + zr * rc.z);
   xr -= (fx * l1.x + fz * l1.z);
   yr -= (fy * l2.y);
   zr -= (fz * l3.z);
   return make_real4(xr, yr, zr, xr * xr + yr * yr + zr * zr);
}


#pragma acc routine seq
__device__
inline real4 image2_orthogonal(real xr, real yr, real zr, real3 l1, real3 l2,
                               real3 l3, real3 ra, real3 rb, real3 rc)
{
   real fx = REAL_FLOOR(0.5f + xr * ra.x);
   real fy = REAL_FLOOR(0.5f + yr * rb.y);
   real fz = REAL_FLOOR(0.5f + zr * rc.z);
   xr -= (fx * l1.x);
   yr -= (fy * l2.y);
   zr -= (fz * l3.z);
   return make_real4(xr, yr, zr, xr * xr + yr * yr + zr * zr);
}


#pragma acc routine seq
template <class T>
__device__
inline real4 image2_triclinic(T r, real3 l1, real3 l2, real3 l3, real3 ra,
                              real3 rb, real3 rc)
{
   return image2_triclinic(r.x, r.y, r.z, l1, l2, l3, ra, rb, rc);
}


#pragma acc routine seq
template <class T>
__device__
inline real4 image2_monoclinic(T r, real3 l1, real3 l2, real3 l3, real3 ra,
                               real3 rb, real3 rc)
{
   return image2_monoclinic(r.x, r.y, r.z, l1, l2, l3, ra, rb, rc);
}


#pragma acc routine seq
template <class T>
__device__
inline real4 image2_orthogonal(T r, real3 l1, real3 l2, real3 l3, real3 ra,
                               real3 rb, real3 rc)
{
   return image2_orthogonal(r.x, r.y, r.z, l1, l2, l3, ra, rb, rc);
}


#ifndef image2
#   define image2(...) image2_triclinic(__VA_ARGS__, TINKER_IMAGE_ARGS)
#endif
TINKER_NAMESPACE_END
