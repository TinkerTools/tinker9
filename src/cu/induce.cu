#include "ff/precision.h"
#include "seq/launch.h"

namespace tinker {
__global__
void pcgUdirV1(int n, const real* restrict polarity, //
   real (*restrict udir)[3], const real (*restrict field)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli = polarity[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         udir[i][j] = poli * field[i][j];
      }
   }
}

__global__
void pcgUdirV2(int n, const real* restrict polarity, real (*restrict udir)[3],
   real (*restrict udirp)[3], const real (*restrict field)[3], const real (*restrict fieldp)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli = polarity[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         udir[i][j] = poli * field[i][j];
         udirp[i][j] = poli * fieldp[i][j];
      }
   }
}
}

namespace tinker {
__global__
void pcgRsd0V1(int n, const real* restrict polarity_inv, real (*restrict rsd)[3], //
   const real (*restrict udir)[3], const real (*restrict uind)[3], const real (*restrict field)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli_inv = polarity_inv[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j)
         rsd[i][j] = (udir[i][j] - uind[i][j]) * poli_inv + field[i][j];
   }
}

__global__
void pcgRsd0V2(int n, const real* restrict polarity_inv, real (*restrict rsd)[3],
   real (*restrict rsp)[3], //
   const real (*restrict udir)[3], const real (*restrict udip)[3], const real (*restrict uind)[3],
   const real (*restrict uinp)[3], const real (*restrict field)[3], const real (*restrict fielp)[3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli_inv = polarity_inv[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         rsd[i][j] = (udir[i][j] - uind[i][j]) * poli_inv + field[i][j];
         rsp[i][j] = (udip[i][j] - uinp[i][j]) * poli_inv + fielp[i][j];
      }
   }
}

__global__
void pcgRsd0V3(int n, const real* restrict polarity_inv, real (*restrict rsd)[3], //
   const real (*restrict udir)[3], const real (*restrict uind)[3], const real (*restrict field)[3],
   const real (*restrict polscale)[3][3])
{
   for (int i = ITHREAD; i < n; i += STRIDE) {
      real poli_inv = polarity_inv[i];
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
         rsd[i][j] = (udir[i][j] - uind[i][0] * polscale[i][0][j] - uind[i][1] * polscale[i][1][j] -
                        uind[i][2] * polscale[i][2][j]) *
               poli_inv +
            field[i][j];
      }
   }
}
}
