#pragma once
#include "add.h"
#include "macro.h"
#include "macro_void_cuda_def.h"


TINKER_NAMESPACE_BEGIN
static constexpr int MAX_BSORDER = 5;


#define bsbuild(j, i) bsbuild_[((i)-1) * bsorder + (j)-1]


// see also subroutine bsplgen in pmestuf.f
#ifdef __CUDACC__
template <int LEVEL, int bsorder>
__device__
void bsplgen(real w, real* restrict thetai, real* restrict bsbuild_)
#else
#pragma acc routine seq
template <int LEVEL>
void bsplgen(real w, real* restrict thetai, int bsorder)
#endif
{
   // e.g. bsorder = 5, theta = T, bsbuild = B

   // LEVEL + 1 <= bsorder

   // LEVEL = 1
   // T(1,1) = B(5,1)
   // T(1,2) = B(5,2)
   // T(1,3) = B(5,3)
   // T(1,4) = B(5,4)
   // T(1,5) = B(5,5)

   // LEVEL = 2
   // T(2,1) = B(4,1)
   // T(2,2) = B(4,2)
   // T(2,3) = B(4,3)
   // T(2,4) = B(4,4)
   // T(2,5) = B(4,5)
   // AND ALL LEVEL = 1

   // LEVEL = 3
   // T(3,1) = B(3,1)
   // T(3,2) = B(3,2)
   // T(3,3) = B(3,3)
   // T(3,4) = B(3,4)
   // T(3,5) = B(3,5)
   // AND ALL LEVEL = 2, 1

   // LEVEL = 4
   // T(4,1) = B(2,1)
   // T(4,2) = B(2,2)
   // T(4,3) = B(2,3)
   // T(4,4) = B(2,4)
   // T(4,5) = B(2,5)
   // AND ALL LEVEL = 3, 2, 1

#ifndef __CUDACC__
   real bsbuild_[MAX_BSORDER * MAX_BSORDER];
#endif

   // initialization to get to 2nd order recursion

   bsbuild(2, 2) = w;
   bsbuild(2, 1) = 1 - w;

   // perform one pass to get to 3rd order recursion

   bsbuild(3, 3) = 0.5f * w * bsbuild(2, 2);
   bsbuild(3, 2) = 0.5f * ((1 + w) * bsbuild(2, 1) + (2 - w) * bsbuild(2, 2));
   bsbuild(3, 1) = 0.5f * (1 - w) * bsbuild(2, 1);

   // compute standard B-spline recursion to desired order

   for (int i = 4; i <= bsorder; ++i) {
      int k = i - 1;
      real denom = REAL_RECIP(k);
      bsbuild(i, i) = denom * w * bsbuild(k, k);
      for (int j = 1; j <= i - 2; j++) {
         bsbuild(i, i - j) = denom *
            ((w + j) * bsbuild(k, i - j - 1) + (i - j - w) * bsbuild(k, i - j));
      }
      bsbuild(i, 1) = denom * (1 - w) * bsbuild(k, 1);
   }

   int k;

   if CONSTEXPR (LEVEL >= 2) {

      // get coefficients for the B-spline first derivative

      k = bsorder - 1;
      bsbuild(k, bsorder) = bsbuild(k, bsorder - 1);
      for (int i = bsorder - 1; i >= 2; --i) {
         bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
      }
      bsbuild(k, 1) = -bsbuild(k, 1);
   }

   if CONSTEXPR (LEVEL >= 3) {

      // get coefficients for the B-spline second derivative

      k = bsorder - 2;
      bsbuild(k, bsorder - 1) = bsbuild(k, bsorder - 2);
      for (int i = bsorder - 2; i >= 2; --i) {
         bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
      }
      bsbuild(k, 1) = -bsbuild(k, 1);
      bsbuild(k, bsorder) = bsbuild(k, bsorder - 1);
      for (int i = bsorder - 1; i >= 2; --i) {
         bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
      }
      bsbuild(k, 1) = -bsbuild(k, 1);
   }

   if CONSTEXPR (LEVEL == 4) {

      // get coefficients for the B-spline third derivative

      k = bsorder - 3;
      bsbuild(k, bsorder - 2) = bsbuild(k, bsorder - 3);
      for (int i = bsorder - 3; i >= 2; --i) {
         bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
      }
      bsbuild(k, 1) = -bsbuild(k, 1);
      bsbuild(k, bsorder - 1) = bsbuild(k, bsorder - 2);
      for (int i = bsorder - 2; i >= 2; --i) {
         bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
      }
      bsbuild(k, 1) = -bsbuild(k, 1);
      bsbuild(k, bsorder) = bsbuild(k, bsorder - 1);
      for (int i = bsorder - 1; i >= 2; --i)
         bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
      bsbuild(k, 1) = -bsbuild(k, 1);
   }

   // copy coefficients from temporary to permanent storage

   for (int i = 1; i <= bsorder; ++i) {
      // keep this line to work on Tesla
      #pragma acc loop seq
      for (int j = 1; j <= LEVEL; ++j) {
         thetai[4 * (i - 1) + (j - 1)] = bsbuild(bsorder - j + 1, i);
      }
   }
}


#undef bsbuild
TINKER_NAMESPACE_END
