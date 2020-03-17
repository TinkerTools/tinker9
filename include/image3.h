#pragma once
#include "box.h"
#include "mathfunc.h"
#include "seq_def.h"


TINKER_NAMESPACE_BEGIN
// SEQ_ROUTINE
// inline void image_triclinic(real& restrict xr, real& restrict yr,
//                             real& restrict zr, real3 l1, real3 l2, real3 l3,
//                             real3 ra, real3 rb, real3 rc)
#define image_triclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc)                    \
   {                                                                           \
      real fx = REAL_FLOOR(0.5f + zr * ra.z + yr * ra.y + xr * ra.x);          \
      real fy = REAL_FLOOR(0.5f + zr * rb.z + yr * rb.y);                      \
      real fz = REAL_FLOOR(0.5f + zr * rc.z);                                  \
      xr -= (fz * l1.z + fy * l1.y + fx * l1.x);                               \
      yr -= (fz * l2.z + fy * l2.y);                                           \
      zr -= (fz * l3.z);                                                       \
   }


// SEQ_ROUTINE
// inline void image_monoclinic(real& restrict xr, real& restrict yr,
//                              real& restrict zr, real3 l1, real3 l2, real3 l3,
//                              real3 ra, real3 rb, real3 rc)
#define image_monoclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc)                   \
   {                                                                           \
      real fx = REAL_FLOOR(0.5f + zr * ra.z + xr * ra.x);                      \
      real fy = REAL_FLOOR(0.5f + yr * rb.y);                                  \
      real fz = REAL_FLOOR(0.5f + zr * rc.z);                                  \
      xr -= (fz * l1.z + fx * l1.x);                                           \
      yr -= (fy * l2.y);                                                       \
      zr -= (fz * l3.z);                                                       \
   }


// SEQ_ROUTINE
// inline void image_orthogonal(real& restrict xr, real& restrict yr,
//                              real& restrict zr, real3 l1, real3 l2, real3 l3,
//                              real3 ra, real3 rb, real3 rc)
#define image_orthogonal(xr, yr, zr, l1, l2, l3, ra, rb, rc)                   \
   {                                                                           \
      real fx = REAL_FLOOR(0.5f + xr * ra.x);                                  \
      real fy = REAL_FLOOR(0.5f + yr * rb.y);                                  \
      real fz = REAL_FLOOR(0.5f + zr * rc.z);                                  \
      xr -= (fx * l1.x);                                                       \
      yr -= (fy * l2.y);                                                       \
      zr -= (fz * l3.z);                                                       \
   }


#define image_oct_macro(xr, yr, zr, l1x, rax)                                  \
   {                                                                           \
      real fx = REAL_FLOOR(0.5f + xr * rax);                                   \
      real fy = REAL_FLOOR(0.5f + yr * rax);                                   \
      real fz = REAL_FLOOR(0.5f + zr * rax);                                   \
      if (REAL_ABS(fx) + REAL_ABS(fy) + REAL_ABS(fz) > 0.75f) {                \
         fx -= REAL_SIGN(0.5, fx);                                             \
         fy -= REAL_SIGN(0.5, fy);                                             \
         fz -= REAL_SIGN(0.5, fz);                                             \
      }                                                                        \
      xr -= (fx * l1x);                                                        \
      yr -= (fy * l1x);                                                        \
      zr -= (fz * l1x);                                                        \
   }
SEQ_ROUTINE
inline void image_oct_devker(real& restrict xr, real& restrict yr,
                             real& restrict zr, real l1x, real rax)
{
   image_oct_macro(xr, yr, zr, l1x, rax);
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
   } else if (l1.x != 0) {
      image_triclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else {
   }
}


SEQ_ROUTINE
inline void image_general(real& restrict xr, real& restrict yr,
                          real& restrict zr, BoxShape sh, real3 l1, real3 l2,
                          real3 l3, real3 ra, real3 rb, real3 rc)
{
   if (sh == ORTHO_BOX) {
      image_orthogonal(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == MONO_BOX) {
      image_monoclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == TRI_BOX) {
      image_triclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == OCT_BOX) {
      image_oct_macro(xr, yr, zr, l1.x, ra.x);
   } else {
   }
}


SEQ_ROUTINE
inline real image2_general(real& restrict xr, real& restrict yr,
                           real& restrict zr, real3 l1, real3 l2, real3 l3,
                           real3 ra, real3 rb, real3 rc)
{
   if (l1.z == 0) {
      image_orthogonal(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (l1.y == 0) {
      image_monoclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (l1.x != 0) {
      image_triclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else {
   }
   return xr * xr + yr * yr + zr * zr;
}


SEQ_ROUTINE
inline real image2_general(real& restrict xr, real& restrict yr,
                           real& restrict zr, BoxShape sh, real3 l1, real3 l2,
                           real3 l3, real3 ra, real3 rb, real3 rc)
{
   if (sh == ORTHO_BOX) {
      image_orthogonal(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == MONO_BOX) {
      image_monoclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == TRI_BOX) {
      image_triclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == OCT_BOX) {
      image_oct_macro(xr, yr, zr, l1.x, ra.x);
   } else {
   }
   return xr * xr + yr * yr + zr * zr;
}


SEQ_ROUTINE
inline real imagen2_general(real& xr, real& yr, real& zr, real3 l1, real3 l2,
                            real3 l3, real3 ra, real3 rb, real3 rc)
{
   if (l1.z == 0) {
      image_orthogonal(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (l1.y == 0) {
      image_monoclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (l1.x != 0) {
      image_triclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else {
   }
   return xr * xr + yr * yr + zr * zr;
}


SEQ_ROUTINE
inline real imagen2_general(real& xr, real& yr, real& zr, BoxShape sh, real3 l1,
                            real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
{
   if (sh == ORTHO_BOX) {
      image_orthogonal(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == MONO_BOX) {
      image_monoclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == TRI_BOX) {
      image_triclinic(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == OCT_BOX) {
      image_oct_macro(xr, yr, zr, l1.x, ra.x);
   } else {
   }
   return xr * xr + yr * yr + zr * zr;
}
TINKER_NAMESPACE_END
