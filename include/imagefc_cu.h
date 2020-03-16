#pragma once
#include "mathfunc.h"


TINKER_NAMESPACE_BEGIN
__device__
inline real3 ftoc_triclinic(real3 f, real3 l1, real3 l2, real3 l3)
{
   f.x = f.z * l1.z + f.y * l1.y + f.x * l1.x;
   f.y = f.z * l2.z + f.y * l2.y;
   f.z = f.z * l3.z;
   return f;
}


__device__
inline real3 ftoc_monoclinic(real3 f, real3 l1, real3 l2, real3 l3)
{
   f.x = f.z * l1.z + f.x * l1.x;
   f.y = f.y * l2.y;
   f.z = f.z * l3.z;
   return f;
}


__device__
inline real3 ftoc_orthogonal(real3 f, real3 l1, real3 l2, real3 l3)
{
   f.x = f.x * l1.x;
   f.y = f.y * l2.y;
   f.z = f.z * l3.z;
   return f;
}


__device__
inline real3 ftoc_general(real3 f, real3 l1, real3 l2, real3 l3)
{
   if (l1.z == 0) {
      return ftoc_orthogonal(f, l1, l2, l3);
   } else if (l1.y == 0) {
      return ftoc_monoclinic(f, l1, l2, l3);
   } else {
      return ftoc_triclinic(f, l1, l2, l3);
   }
}


#ifndef ftoc
#   define ftoc(f) ftoc_general(f, lvec1, lvec2, lvec3)
#endif


//====================================================================//


__device__
inline real3 imagef(real3 f)
{
   f.x -= REAL_FLOOR(0.5f + f.x);
   f.y -= REAL_FLOOR(0.5f + f.y);
   f.z -= REAL_FLOOR(0.5f + f.z);
   return f;
}


//====================================================================//


__device__
inline real3 ctof_triclinic(real xr, real yr, real zr, real3 ra, real3 rb,
                            real3 rc)
{
   real3 f;
   f.x = zr * ra.z + yr * ra.y + xr * ra.x;
   f.y = zr * rb.z + yr * rb.y;
   f.z = zr * rc.z;
   return f;
}


__device__
inline real3 ctof_monoclinic(real xr, real yr, real zr, real3 ra, real3 rb,
                             real3 rc)
{
   real3 f;
   f.x = zr * ra.z + xr * ra.x;
   f.y = yr * rb.y;
   f.z = zr * rc.z;
   return f;
}


__device__
inline real3 ctof_orthogonal(real xr, real yr, real zr, real3 ra, real3 rb,
                             real3 rc)
{
   real3 f;
   f.x = xr * ra.x;
   f.y = yr * rb.y;
   f.z = zr * rc.z;
   return f;
}


__device__
inline real3 imagectof_general(real xr, real yr, real zr, real3 ra, real3 rb,
                               real3 rc)
{
   if (ra.z == 0) {
      return imagef(ctof_orthogonal(xr, yr, zr, ra, rb, rc));
   } else if (ra.y == 0) {
      return imagef(ctof_monoclinic(xr, yr, zr, ra, rb, rc));
   } else {
      return imagef(ctof_triclinic(xr, yr, zr, ra, rb, rc));
   }
}


#ifndef imagectof
#   define imagectof(xr, yr, zr)                                               \
      imagectof_general(xr, yr, zr, recipa, recipb, recipc)
#endif
TINKER_NAMESPACE_END
