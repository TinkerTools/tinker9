#pragma once
#include "macro.h"


TINKER_NAMESPACE_BEGIN
__global__
static void epolar_trq_cu(real* restrict trqx, real* restrict trqy,
                          real* restrict trqz, int n,
                          const real (*restrict rpole)[10],
                          const real (*restrict ufld)[3],
                          const real (*restrict dufld)[6])
{
   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      real dix = rpole[i][mpl_pme_x];
      real diy = rpole[i][mpl_pme_y];
      real diz = rpole[i][mpl_pme_z];
      real qixx = rpole[i][mpl_pme_xx];
      real qixy = rpole[i][mpl_pme_xy];
      real qixz = rpole[i][mpl_pme_xz];
      real qiyy = rpole[i][mpl_pme_yy];
      real qiyz = rpole[i][mpl_pme_yz];
      real qizz = rpole[i][mpl_pme_zz];


      real tep1 = diz * ufld[i][1] - diy * ufld[i][2] + qixz * dufld[i][1] -
         qixy * dufld[i][3] + 2 * qiyz * (dufld[i][2] - dufld[i][5]) +
         (qizz - qiyy) * dufld[i][4];
      real tep2 = dix * ufld[i][2] - diz * ufld[i][0] - qiyz * dufld[i][1] +
         qixy * dufld[i][4] + 2 * qixz * (dufld[i][5] - dufld[i][0]) +
         (qixx - qizz) * dufld[i][3];
      real tep3 = diy * ufld[i][0] - dix * ufld[i][1] + qiyz * dufld[i][3] -
         qixz * dufld[i][4] + 2 * qixy * (dufld[i][0] - dufld[i][2]) +
         (qiyy - qixx) * dufld[i][1];


      trqx[i] += tep1;
      trqy[i] += tep2;
      trqz[i] += tep3;
   }
}
TINKER_NAMESPACE_END
