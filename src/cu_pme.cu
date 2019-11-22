#include "add.cuh"
#include "box.h"
#include "e_mpole.h"
#include "error.h"
#include "launch.cuh"
#include "md.h"
#include "pme.h"
#include "seq_pme.h"


TINKER_NAMESPACE_BEGIN
static constexpr int PME_BLOCKDIM = 64;


enum
{
   PCHG_GRID = 1,
   MPOLE_GRID,
   UIND_GRID,
   UIND_GRID_FPHI2,
   DISP_GRID
};


template <int WHAT, int bsorder>
__global__
void grid_tmpl_cu(const real* restrict x, const real* restrict y,
                  const real* restrict z, int n, int nfft1, int nfft2,
                  int nfft3, const real (*fmp)[10], real* restrict qgrid,
                  const Box* box)
{
   real thetai1[4 * MAX_BSORDER];
   real thetai2[4 * MAX_BSORDER];
   real thetai3[4 * MAX_BSORDER];
   __shared__ real sharedarray[MAX_BSORDER * MAX_BSORDER * PME_BLOCKDIM];
   real* restrict array = &sharedarray[MAX_BSORDER * MAX_BSORDER * threadIdx.x];


   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      real w1 =
         xi * box->recip[0][0] + yi * box->recip[0][1] + zi * box->recip[0][2];
      w1 = w1 + 0.5f - REAL_FLOOR(w1 + 0.5f);
      real fr1 = nfft1 * w1;
      int igrid1 = REAL_FLOOR(fr1);
      w1 = fr1 - igrid1;


      real w2 =
         xi * box->recip[1][0] + yi * box->recip[1][1] + zi * box->recip[1][2];
      w2 = w2 + 0.5f - REAL_FLOOR(w2 + 0.5f);
      real fr2 = nfft2 * w2;
      int igrid2 = REAL_FLOOR(fr2);
      w2 = fr2 - igrid2;


      real w3 =
         xi * box->recip[2][0] + yi * box->recip[2][1] + zi * box->recip[2][2];
      w3 = w3 + 0.5f - REAL_FLOOR(w3 + 0.5f);
      real fr3 = nfft3 * w3;
      int igrid3 = REAL_FLOOR(fr3);
      w3 = fr3 - igrid3;


      igrid1 = igrid1 - bsorder + 1;
      igrid2 = igrid2 - bsorder + 1;
      igrid3 = igrid3 - bsorder + 1;
      igrid1 += (igrid1 < 0 ? nfft1 : 0);
      igrid2 += (igrid2 < 0 ? nfft2 : 0);
      igrid3 += (igrid3 < 0 ? nfft3 : 0);


      if CONSTEXPR (WHAT == MPOLE_GRID) {
         bsplgen<3, bsorder>(w1, thetai1, array);
         bsplgen<3, bsorder>(w2, thetai2, array);
         bsplgen<3, bsorder>(w3, thetai3, array);
      }


      if CONSTEXPR (WHAT == MPOLE_GRID) {
         real fmpi0 = fmp[i][mpl_pme_0];
         real fmpix = fmp[i][mpl_pme_x];
         real fmpiy = fmp[i][mpl_pme_y];
         real fmpiz = fmp[i][mpl_pme_z];
         real fmpixx = fmp[i][mpl_pme_xx];
         real fmpiyy = fmp[i][mpl_pme_yy];
         real fmpizz = fmp[i][mpl_pme_zz];
         real fmpixy = fmp[i][mpl_pme_xy];
         real fmpixz = fmp[i][mpl_pme_xz];
         real fmpiyz = fmp[i][mpl_pme_yz];
         for (int iz = 0; iz < bsorder; ++iz) {
            int zbase = igrid3 + iz;
            zbase -= (zbase >= nfft3 ? nfft3 : 0);
            zbase *= (nfft1 * nfft2);
            real v0 = thetai3[4 * iz];
            real v1 = thetai3[4 * iz + 1];
            real v2 = thetai3[4 * iz + 2];
            for (int iy = 0; iy < bsorder; ++iy) {
               int ybase = igrid2 + iy;
               ybase -= (ybase >= nfft2 ? nfft2 : 0);
               ybase *= nfft1;
               real u0 = thetai2[4 * iy];
               real u1 = thetai2[4 * iy + 1];
               real u2 = thetai2[4 * iy + 2];
               // fmp: 0, x, y, z, xx, yy, zz, xy, xz, yz
               //      1, 2, 3, 4,  5,  6,  7,  8,  9, 10
               real term0 = fmpi0 * u0 * v0 + fmpiy * u1 * v0 +
                  fmpiz * u0 * v1 + fmpiyy * u2 * v0 + fmpizz * u0 * v2 +
                  fmpiyz * u1 * v1;
               real term1 =
                  fmpix * u0 * v0 + fmpixy * u1 * v0 + fmpixz * u0 * v1;
               real term2 = fmpixx * u0 * v0;
               for (int ix = 0; ix < bsorder; ++ix) {
                  int xbase = igrid1 + ix;
                  xbase -= (xbase >= nfft1 ? nfft1 : 0);
                  int index = xbase + ybase + zbase;
                  real t0 = thetai1[4 * ix];
                  real t1 = thetai1[4 * ix + 1];
                  real t2 = thetai1[4 * ix + 2];
                  atomic_add(term0 * t0 + term1 * t1 + term2 * t2, qgrid,
                             2 * index);
               }
            } // end for (int iy)
         }
      } // end if (WHAT == MPOLE_GRID)
   }
}


void grid_mpole_cu(PMEUnit pme_u, real (*fmp)[10])
{
   auto& st = *pme_u;
   int n1 = st.nfft1;
   int n2 = st.nfft2;
   int n3 = st.nfft3;
   int nt = n1 * n2 * n3;


   if (st.bsorder != 5)
      TINKER_THROW(
         format("grid_mpole_cu(): bsorder is {}; must be 5.\n", st.bsorder));


   device_array::zero(2 * nt, st.qgrid);
   auto ker = grid_tmpl_cu<MPOLE_GRID, 5>;
   launch_kernel2(PME_BLOCKDIM, n, ker, x, y, z, n, n1, n2, n3, fmp, st.qgrid,
                  box);
}
TINKER_NAMESPACE_END
