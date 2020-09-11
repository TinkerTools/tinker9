#pragma once
#include "box.h"
#include "mathfunc.h"
#include "seq_def.h"


namespace tinker {
#define IMAGE_TRI__(xr, yr, zr, l1, l2, l3, ra, rb, rc)                        \
   {                                                                           \
      real fx = REAL_FLOOR(0.5f + zr * ra.z + yr * ra.y + xr * ra.x);          \
      real fy = REAL_FLOOR(0.5f + zr * rb.z + yr * rb.y);                      \
      real fz = REAL_FLOOR(0.5f + zr * rc.z);                                  \
      xr -= (fz * l1.z + fy * l1.y + fx * l1.x);                               \
      yr -= (fz * l2.z + fy * l2.y);                                           \
      zr -= (fz * l3.z);                                                       \
   }
#define IMAGE_MONO__(xr, yr, zr, l1, l2, l3, ra, rb, rc)                       \
   {                                                                           \
      real fx = REAL_FLOOR(0.5f + zr * ra.z + xr * ra.x);                      \
      real fy = REAL_FLOOR(0.5f + yr * rb.y);                                  \
      real fz = REAL_FLOOR(0.5f + zr * rc.z);                                  \
      xr -= (fz * l1.z + fx * l1.x);                                           \
      yr -= (fy * l2.y);                                                       \
      zr -= (fz * l3.z);                                                       \
   }
#define IMAGE_ORTHO__(xr, yr, zr, l1, l2, l3, ra, rb, rc)                      \
   {                                                                           \
      real fx = REAL_FLOOR(0.5f + xr * ra.x);                                  \
      real fy = REAL_FLOOR(0.5f + yr * rb.y);                                  \
      real fz = REAL_FLOOR(0.5f + zr * rc.z);                                  \
      xr -= (fx * l1.x);                                                       \
      yr -= (fy * l2.y);                                                       \
      zr -= (fz * l3.z);                                                       \
   }
#define IMAGE_OCT__(xr, yr, zr, l1, l2, l3, ra, rb, rc)                        \
   {                                                                           \
      real fx = xr * ra.x;                                                     \
      real fy = yr * ra.x;                                                     \
      real fz = zr * ra.x;                                                     \
      fx -= REAL_FLOOR(0.5f + fx);                                             \
      fy -= REAL_FLOOR(0.5f + fy);                                             \
      fz -= REAL_FLOOR(0.5f + fz);                                             \
      if (REAL_ABS(fx) + REAL_ABS(fy) + REAL_ABS(fz) > 0.75f) {                \
         fx -= REAL_SIGN(0.5f, fx);                                            \
         fy -= REAL_SIGN(0.5f, fy);                                            \
         fz -= REAL_SIGN(0.5f, fz);                                            \
      }                                                                        \
      xr = fx * l1.x;                                                          \
      yr = fy * l1.x;                                                          \
      zr = fz * l1.x;                                                          \
   }


struct PBC_ORTHO
{
   SEQ_ROUTINE
   static void img(real& restrict xr, real& restrict yr, real& restrict zr,
                   real3 l1, real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
   {
      IMAGE_ORTHO__(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   }


   SEQ_ROUTINE
   static real img2(real& restrict xr, real& restrict yr, real& restrict zr,
                    real3 l1, real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
   {
      img(xr, yr, zr, l1, l2, l3, ra, rb, rc);
      return xr * xr + yr * yr + zr * zr;
   }
};


struct PBC_MONO
{
   SEQ_ROUTINE
   static void img(real& restrict xr, real& restrict yr, real& restrict zr,
                   real3 l1, real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
   {
      IMAGE_MONO__(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   }


   SEQ_ROUTINE
   static real img2(real& restrict xr, real& restrict yr, real& restrict zr,
                    real3 l1, real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
   {
      img(xr, yr, zr, l1, l2, l3, ra, rb, rc);
      return xr * xr + yr * yr + zr * zr;
   }
};


struct PBC_TRI
{
   SEQ_ROUTINE
   static void img(real& restrict xr, real& restrict yr, real& restrict zr,
                   real3 l1, real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
   {
      IMAGE_TRI__(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   }


   SEQ_ROUTINE
   static real img2(real& restrict xr, real& restrict yr, real& restrict zr,
                    real3 l1, real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
   {
      img(xr, yr, zr, l1, l2, l3, ra, rb, rc);
      return xr * xr + yr * yr + zr * zr;
   }
};


struct PBC_OCT
{
   SEQ_ROUTINE
   static void img(real& restrict xr, real& restrict yr, real& restrict zr,
                   real3 l1, real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
   {
      IMAGE_OCT__(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   }


   SEQ_ROUTINE
   static real img2(real& restrict xr, real& restrict yr, real& restrict zr,
                    real3 l1, real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
   {
      img(xr, yr, zr, l1, l2, l3, ra, rb, rc);
      return xr * xr + yr * yr + zr * zr;
   }
};


//====================================================================//


SEQ_ROUTINE
inline void image_general(real& restrict xr, real& restrict yr,
                          real& restrict zr, BoxShape sh, real3 l1, real3 l2,
                          real3 l3, real3 ra, real3 rb, real3 rc)
{
   if (sh == ORTHO_BOX) {
      IMAGE_ORTHO__(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == MONO_BOX) {
      IMAGE_MONO__(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == TRI_BOX) {
      IMAGE_TRI__(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == OCT_BOX) {
      IMAGE_OCT__(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else {
      // UNBOND_BOX
   }
}


SEQ_ROUTINE
inline real image2_general(real& restrict xr, real& restrict yr,
                           real& restrict zr, BoxShape sh, real3 l1, real3 l2,
                           real3 l3, real3 ra, real3 rb, real3 rc)
{
   image_general(xr, yr, zr, sh, l1, l2, l3, ra, rb, rc);
   return xr * xr + yr * yr + zr * zr;
}


SEQ_ROUTINE
inline real imagen2_general(real& xr, real& yr, real& zr, BoxShape sh, real3 l1,
                            real3 l2, real3 l3, real3 ra, real3 rb, real3 rc)
{
   if (sh == ORTHO_BOX) {
      IMAGE_ORTHO__(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == MONO_BOX) {
      IMAGE_MONO__(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == TRI_BOX) {
      IMAGE_TRI__(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else if (sh == OCT_BOX) {
      IMAGE_OCT__(xr, yr, zr, l1, l2, l3, ra, rb, rc);
   } else {
      // UNBOND_BOX
   }
   return xr * xr + yr * yr + zr * zr;
}
}
