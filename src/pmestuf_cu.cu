#include "add.h"
#include "box.h"
#include "empole.h"
#include "launch.h"
#include "md.h"
#include "named_struct.h"
#include "pmestuf.h"
#include "seq_bsplgen.h"
#include "spatial.h"


TINKER_NAMESPACE_BEGIN
// compute theta values on the fly
template <class T, int bsorder>
__global__
void grid_put_cu1(const real* restrict x, const real* restrict y,
                  const real* restrict z, int n, int nfft1, int nfft2,
                  int nfft3, const real* restrict ptr1,
                  const real* restrict ptr2, real* restrict qgrid,
                  real3 recip_a, real3 recip_b, real3 recip_c)
{
   real thetai1[4 * 5];
   real thetai2[4 * 5];
   real thetai3[4 * 5];
   __shared__ real sharedarray[5 * 5 * PME_BLOCKDIM];
   real* restrict array = &sharedarray[5 * 5 * threadIdx.x];


   MAYBE_UNUSED const real* pchg = ptr1;
   MAYBE_UNUSED const real(*fmp)[10] = (real(*)[10])ptr1;
   MAYBE_UNUSED const real(*fuind)[3] = (real(*)[3])ptr1;
   MAYBE_UNUSED const real(*fuinp)[3] = (real(*)[3])ptr2;


   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      real w1 = xi * recip_a.x + yi * recip_a.y + zi * recip_a.z;
      w1 = w1 + 0.5f - REAL_FLOOR(w1 + 0.5f);
      real fr1 = nfft1 * w1;
      int igrid1 = REAL_FLOOR(fr1);
      w1 = fr1 - igrid1;


      real w2 = xi * recip_b.x + yi * recip_b.y + zi * recip_b.z;
      w2 = w2 + 0.5f - REAL_FLOOR(w2 + 0.5f);
      real fr2 = nfft2 * w2;
      int igrid2 = REAL_FLOOR(fr2);
      w2 = fr2 - igrid2;


      real w3 = xi * recip_c.x + yi * recip_c.y + zi * recip_c.z;
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


      if CONSTEXPR (eq<T, PCHG>()) {
         real chgi = pchg[i];
         if (chgi == 0)
            continue;


         bsplgen<1, bsorder>(w1, thetai1, array);
         bsplgen<1, bsorder>(w2, thetai2, array);
         bsplgen<1, bsorder>(w3, thetai3, array);


         for (int iz = 0; iz < bsorder; ++iz) {
            int zbase = igrid3 + iz;
            zbase -= (zbase >= nfft3 ? nfft3 : 0);
            zbase *= (nfft1 * nfft2);
            real v0 = thetai3[4 * iz] * chgi;
            for (int iy = 0; iy < bsorder; ++iy) {
               int ybase = igrid2 + iy;
               ybase -= (ybase >= nfft2 ? nfft2 : 0);
               ybase *= nfft1;
               real u0 = thetai2[4 * iy] * v0;
               for (int ix = 0; ix < bsorder; ++ix) {
                  int xbase = igrid1 + ix;
                  xbase -= (xbase >= nfft1 ? nfft1 : 0);
                  int index = xbase + ybase + zbase;
                  real term = thetai1[4 * ix] * u0;
                  atomic_add(term, qgrid, 2 * index);
               }
            }
         }
      } // end if (PCHG)


      if CONSTEXPR (eq<T, MPOLE>()) {
         bsplgen<3, bsorder>(w1, thetai1, array);
         bsplgen<3, bsorder>(w2, thetai2, array);
         bsplgen<3, bsorder>(w3, thetai3, array);


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
      } // end if (MPOLE)


      if CONSTEXPR (eq<T, UIND>()) {
         bsplgen<2, bsorder>(w1, thetai1, array);
         bsplgen<2, bsorder>(w2, thetai2, array);
         bsplgen<2, bsorder>(w3, thetai3, array);


         real fuindi0 = fuind[i][0];
         real fuindi1 = fuind[i][1];
         real fuindi2 = fuind[i][2];
         real fuinpi0 = fuinp[i][0];
         real fuinpi1 = fuinp[i][1];
         real fuinpi2 = fuinp[i][2];
         #pragma unroll
         for (int iz = 0; iz < bsorder; ++iz) {
            int zbase = igrid3 + iz;
            zbase -= (zbase >= nfft3 ? nfft3 : 0);
            zbase *= (nfft1 * nfft2);
            real v0 = thetai3[4 * iz];
            real v1 = thetai3[4 * iz + 1];
            #pragma unroll
            for (int iy = 0; iy < bsorder; ++iy) {
               int ybase = igrid2 + iy;
               ybase -= (ybase >= nfft2 ? nfft2 : 0);
               ybase *= nfft1;
               real u0 = thetai2[4 * iy];
               real u1 = thetai2[4 * iy + 1];
               real term01 = fuindi1 * u1 * v0 + fuindi2 * u0 * v1;
               real term11 = fuindi0 * u0 * v0;
               real term02 = fuinpi1 * u1 * v0 + fuinpi2 * u0 * v1;
               real term12 = fuinpi0 * u0 * v0;
               #pragma unroll
               for (int ix = 0; ix < bsorder; ++ix) {
                  int xbase = igrid1 + ix;
                  xbase -= (xbase >= nfft1 ? nfft1 : 0);
                  int index = xbase + ybase + zbase;
                  real t0 = thetai1[4 * ix];
                  real t1 = thetai1[4 * ix + 1];
                  atomic_add(term01 * t0 + term11 * t1, qgrid, 2 * index);
                  atomic_add(term02 * t0 + term12 * t1, qgrid, 2 * index + 1);
               }
            } // end for (int iy)
         }
      } // end if (UIND)
   }
}


// use pre-computed theta values
template <class T, int bsorder>
__global__
void grid_put_cu2(const int* restrict igrid, const real* restrict thetai1,
                  const real* restrict thetai2, const real* restrict thetai3,
                  const Spatial::SortedAtom* restrict sorted, int n,
                  int padded_n, int nfft1, int nfft2, int nfft3,
                  const real* restrict ptr1, const real* ptr2,
                  real* restrict qgrid)
{
   constexpr int bso2 = bsorder * bsorder;
   constexpr int bso3 = bsorder * bso2;
   for (int m = threadIdx.x + blockIdx.x * blockDim.x; m < n * bso3;
        m += blockDim.x * gridDim.x) {
      // m = i0 * bso3 + j;
      int i0 = m / bso3;
      int j = m - i0 * bso3;


      int i = sorted[i0].unsorted;
      int igrid1 = igrid[3 * i + 0];
      int igrid2 = igrid[3 * i + 1];
      int igrid3 = igrid[3 * i + 2];


      int iz = j / bso2;
      j -= iz * bso2;
      int iy = j / bsorder;
      int ix = j - (j / bsorder) * bsorder;


      real v0 = thetai3[(4 * iz + 0) * padded_n + i];
      real v1 = thetai3[(4 * iz + 1) * padded_n + i];
      int zbase = igrid3 + iz;
      zbase -= (zbase >= nfft3 ? nfft3 : 0);
      zbase *= (nfft1 * nfft2);


      real u0 = thetai2[(4 * iy + 0) * padded_n + i];
      real u1 = thetai2[(4 * iy + 1) * padded_n + i];
      int ybase = igrid2 + iy;
      ybase -= (ybase >= nfft2 ? nfft2 : 0);
      ybase *= nfft1;


      real t0 = thetai1[(4 * ix + 0) * padded_n + i];
      real t1 = thetai1[(4 * ix + 1) * padded_n + i];
      int xbase = igrid1 + ix;
      xbase -= (xbase >= nfft1 ? nfft1 : 0);
      int index = xbase + ybase + zbase;


      if CONSTEXPR (eq<T, MPOLE>()) {
         real v2 = thetai3[(4 * iz + 2) * padded_n + i];
         real u2 = thetai2[(4 * iy + 2) * padded_n + i];
         real t2 = thetai1[(4 * ix + 2) * padded_n + i];
         real fmpi0 = ptr1[i * 10 + mpl_pme_0];
         real fmpix = ptr1[i * 10 + mpl_pme_x];
         real fmpiy = ptr1[i * 10 + mpl_pme_y];
         real fmpiz = ptr1[i * 10 + mpl_pme_z];
         real fmpixx = ptr1[i * 10 + mpl_pme_xx];
         real fmpiyy = ptr1[i * 10 + mpl_pme_yy];
         real fmpizz = ptr1[i * 10 + mpl_pme_zz];
         real fmpixy = ptr1[i * 10 + mpl_pme_xy];
         real fmpixz = ptr1[i * 10 + mpl_pme_xz];
         real fmpiyz = ptr1[i * 10 + mpl_pme_yz];
         real term0 = fmpi0 * u0 * v0 + fmpiy * u1 * v0 + fmpiz * u0 * v1 +
            fmpiyy * u2 * v0 + fmpizz * u0 * v2 + fmpiyz * u1 * v1;
         real term1 = fmpix * u0 * v0 + fmpixy * u1 * v0 + fmpixz * u0 * v1;
         real term2 = fmpixx * u0 * v0;
         atomic_add(term0 * t0 + term1 * t1 + term2 * t2, qgrid, 2 * index);
      }


      if CONSTEXPR (eq<T, UIND>()) {
         real3 fd =
            make_real3(ptr1[3 * i + 0], ptr1[3 * i + 1], ptr1[3 * i + 2]);
         real3 fp =
            make_real3(ptr2[3 * i + 0], ptr2[3 * i + 1], ptr2[3 * i + 2]);
         real3 tuv = make_real3(t1 * u0 * v0, t0 * u1 * v0, t0 * u0 * v1);
         atomic_add(dot3(fd, tuv), qgrid, 2 * index);
         atomic_add(dot3(fp, tuv), qgrid, 2 * index + 1);
      }
   }
}


void grid_pchg_cu(PMEUnit pme_u, real* pchg)
{
   auto& st = *pme_u;
   int n1 = st.nfft1;
   int n2 = st.nfft2;
   int n3 = st.nfft3;
   int nt = n1 * n2 * n3;


   darray::zero(PROCEED_NEW_Q, 2 * nt, st.qgrid);
   auto ker = grid_put_cu1<PCHG, 5>;
   launch_k2s(nonblk, PME_BLOCKDIM, n, ker, x, y, z, n, n1, n2, n3, pchg,
              nullptr, st.qgrid, recipa, recipb, recipc);
}


void grid_mpole_cu(PMEUnit pme_u, real (*fmp)[10])
{
   auto& st = *pme_u;
   int n1 = st.nfft1;
   int n2 = st.nfft2;
   int n3 = st.nfft3;
   int nt = n1 * n2 * n3;


   darray::zero(PROCEED_NEW_Q, 2 * nt, st.qgrid);
   if (TINKER_CU_THETA_ON_THE_FLY_GRID_MPOLE) {
      auto ker = grid_put_cu1<MPOLE, 5>;
      launch_k2s(nonblk, PME_BLOCKDIM, n, ker, x, y, z, n, n1, n2, n3,
                 (const real*)fmp, nullptr, st.qgrid, recipa, recipb, recipc);
   } else {
      auto ker = grid_put_cu2<MPOLE, 5>;
      int npa = 5 * 5 * 5 * n;
      launch_k1s(nonblk, npa, ker, st.igrid, st.thetai1, st.thetai2, st.thetai3,
                 mspatial_unit->sorted, n, padded_n, n1, n2, n3,
                 (const real*)fmp, nullptr, st.qgrid);
   }
}


void grid_uind_cu(PMEUnit pme_u, real (*fuind)[3], real (*fuinp)[3])
{
   auto& st = *pme_u;
   int n1 = st.nfft1;
   int n2 = st.nfft2;
   int n3 = st.nfft3;
   int nt = n1 * n2 * n3;


   darray::zero(PROCEED_NEW_Q, 2 * nt, st.qgrid);
   if (TINKER_CU_THETA_ON_THE_FLY_GRID_UIND) {
      auto ker = grid_put_cu1<UIND, 5>;
      launch_k2s(nonblk, PME_BLOCKDIM, n, ker, x, y, z, n, n1, n2, n3,
                 (const real*)fuind, (const real*)fuinp, st.qgrid, recipa,
                 recipb, recipc);
   } else {
      auto ker = grid_put_cu2<UIND, 5>;
      int npa = 5 * 5 * 5 * n;
      launch_k1s(nonblk, npa, ker, st.igrid, st.thetai1, st.thetai2, st.thetai3,
                 mspatial_unit->sorted, n, padded_n, n1, n2, n3,
                 (const real*)fuind, (const real*)fuinp, st.qgrid);
   }
}


template <int LEVEL, int bsorder>
__global__
void bspline_fill_cu1(int* restrict igrid, real* restrict thetai1,
                      real* restrict thetai2, real* restrict thetai3,
                      const real* restrict x, const real* restrict y,
                      const real* restrict z, int n, int padded_n, int nfft1,
                      int nfft2, int nfft3, real3 recip_a, real3 recip_b,
                      real3 recip_c)
{
   const int nfft4[3] = {nfft1, nfft2, nfft3};
   const real3 recip4[3] = {recip_a, recip_b, recip_c};
   real* const thetai[3] = {thetai1, thetai2, thetai3};
   real array[5 * 5];


   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      int igridi[3];
      for (int j = 0; j < 3; ++j) {
         real w4 = xi * recip4[j].x + yi * recip4[j].y + zi * recip4[j].z;
         w4 = w4 + 0.5f - REAL_FLOOR(w4 + 0.5f);
         real fr4 = nfft4[j] * w4;
         int igrid4 = REAL_FLOOR(fr4);
         w4 = fr4 - igrid4;
         igrid4 = igrid4 - bsorder + 1;
         igrid4 += (igrid4 < 0 ? nfft4[j] : 0);
         // write output
         igridi[j] = igrid4;
         bsplgen2<LEVEL, bsorder>(w4, thetai[j], i, padded_n, array);
      }
      igrid[3 * i + 0] = igridi[0];
      igrid[3 * i + 1] = igridi[1];
      igrid[3 * i + 2] = igridi[2];
   }
}


void bspline_fill_cu(PMEUnit u, int level)
{
   auto& st = *u;
   if (level == 2) {
      auto ker = bspline_fill_cu1<2, 5>;
      launch_k1s(nonblk, n, ker, st.igrid, st.thetai1, st.thetai2, st.thetai3,
                 x, y, z, n, padded_n, st.nfft1, st.nfft2, st.nfft3, recipa,
                 recipb, recipc);
   } else if (level == 3) {
      auto ker = bspline_fill_cu1<3, 5>;
      launch_k1s(nonblk, n, ker, st.igrid, st.thetai1, st.thetai2, st.thetai3,
                 x, y, z, n, padded_n, st.nfft1, st.nfft2, st.nfft3, recipa,
                 recipb, recipc);
   }
}


// compute theta values on the fly
template <class T, int bsorder>
__global__
void fphi_get_cu(int n, int nfft1, int nfft2, int nfft3, const real* restrict x,
                 const real* restrict y, const real* restrict z,
                 real* restrict opt1, real* restrict opt2, real* restrict opt3,
                 const real* restrict qgrid, real3 recip_a, real3 recip_b,
                 real3 recip_c)
{
   real thetai1[4 * 5];
   real thetai2[4 * 5];
   real thetai3[4 * 5];
   real array[5 * 5];


   MAYBE_UNUSED real(*restrict fphi)[20] = (real(*)[20])opt1;
   MAYBE_UNUSED real(*restrict fdip_phi1)[10] = (real(*)[10])opt1;
   MAYBE_UNUSED real(*restrict fdip_phi2)[10] = (real(*)[10])opt2;
   MAYBE_UNUSED real(*restrict fdip_sum_phi)[20] = (real(*)[20])opt3;


   for (int i = threadIdx.x + blockIdx.x * blockDim.x; i < n;
        i += blockDim.x * gridDim.x) {
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      real w1 = xi * recip_a.x + yi * recip_a.y + zi * recip_a.z;
      w1 = w1 + 0.5f - REAL_FLOOR(w1 + 0.5f);
      real fr1 = nfft1 * w1;
      int igrid1 = REAL_FLOOR(fr1);
      w1 = fr1 - igrid1;


      real w2 = xi * xi * recip_b.x + yi * recip_b.y + zi * recip_b.z;
      w2 = w2 + 0.5f - REAL_FLOOR(w2 + 0.5f);
      real fr2 = nfft2 * w2;
      int igrid2 = REAL_FLOOR(fr2);
      w2 = fr2 - igrid2;


      real w3 = xi * recip_c.x + yi * recip_c.y + zi * recip_c.z;
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

      if CONSTEXPR (eq<T, MPOLE>() || eq<T, UIND>() || eq<T, UIND2>()) {
         bsplgen<4, bsorder>(w1, thetai1, array);
         bsplgen<4, bsorder>(w2, thetai2, array);
         bsplgen<4, bsorder>(w3, thetai3, array);
      }


      if CONSTEXPR (eq<T, MPOLE>()) {
         real tuv000 = 0;
         real tuv001 = 0;
         real tuv010 = 0;
         real tuv100 = 0;
         real tuv200 = 0;
         real tuv020 = 0;
         real tuv002 = 0;
         real tuv110 = 0;
         real tuv101 = 0;
         real tuv011 = 0;
         real tuv300 = 0;
         real tuv030 = 0;
         real tuv003 = 0;
         real tuv210 = 0;
         real tuv201 = 0;
         real tuv120 = 0;
         real tuv021 = 0;
         real tuv102 = 0;
         real tuv012 = 0;
         real tuv111 = 0;
         for (int iz = 0; iz < bsorder; ++iz) {
            int zbase = igrid3 + iz;
            zbase -= (zbase >= nfft3 ? nfft3 : 0);
            zbase *= (nfft1 * nfft2);
            real v0 = thetai3[4 * iz];
            real v1 = thetai3[4 * iz + 1];
            real v2 = thetai3[4 * iz + 2];
            real v3 = thetai3[4 * iz + 3];
            real tu00 = 0;
            real tu10 = 0;
            real tu01 = 0;
            real tu20 = 0;
            real tu11 = 0;
            real tu02 = 0;
            real tu30 = 0;
            real tu21 = 0;
            real tu12 = 0;
            real tu03 = 0;
            for (int iy = 0; iy < bsorder; ++iy) {
               int ybase = igrid2 + iy;
               ybase -= (ybase >= nfft2 ? nfft2 : 0);
               ybase *= nfft1;
               real u0 = thetai2[4 * iy];
               real u1 = thetai2[4 * iy + 1];
               real u2 = thetai2[4 * iy + 2];
               real u3 = thetai2[4 * iy + 3];
               real t0 = 0;
               real t1 = 0;
               real t2 = 0;
               real t3 = 0;
               for (int ix = 0; ix < bsorder; ++ix) {
                  int xbase = igrid1 + ix;
                  xbase -= (xbase >= nfft1 ? nfft1 : 0);
                  real tq = qgrid[2 * (xbase + ybase + zbase)];
                  t0 += tq * thetai1[4 * ix];
                  t1 += tq * thetai1[4 * ix + 1];
                  t2 += tq * thetai1[4 * ix + 2];
                  t3 += tq * thetai1[4 * ix + 3];
               }
               tu00 += t0 * u0;
               tu10 += t1 * u0;
               tu01 += t0 * u1;
               tu20 += t2 * u0;
               tu11 += t1 * u1;
               tu02 += t0 * u2;
               tu30 += t3 * u0;
               tu21 += t2 * u1;
               tu12 += t1 * u2;
               tu03 += t0 * u3;
            }
            tuv000 += tu00 * v0;
            tuv100 += tu10 * v0;
            tuv010 += tu01 * v0;
            tuv001 += tu00 * v1;
            tuv200 += tu20 * v0;
            tuv020 += tu02 * v0;
            tuv002 += tu00 * v2;
            tuv110 += tu11 * v0;
            tuv101 += tu10 * v1;
            tuv011 += tu01 * v1;
            tuv300 += tu30 * v0;
            tuv030 += tu03 * v0;
            tuv003 += tu00 * v3;
            tuv210 += tu21 * v0;
            tuv201 += tu20 * v1;
            tuv120 += tu12 * v0;
            tuv021 += tu02 * v1;
            tuv102 += tu10 * v2;
            tuv012 += tu01 * v2;
            tuv111 += tu11 * v1;
         }
         fphi[i][0] = tuv000;
         fphi[i][1] = tuv100;
         fphi[i][2] = tuv010;
         fphi[i][3] = tuv001;
         fphi[i][4] = tuv200;
         fphi[i][5] = tuv020;
         fphi[i][6] = tuv002;
         fphi[i][7] = tuv110;
         fphi[i][8] = tuv101;
         fphi[i][9] = tuv011;
         fphi[i][10] = tuv300;
         fphi[i][11] = tuv030;
         fphi[i][12] = tuv003;
         fphi[i][13] = tuv210;
         fphi[i][14] = tuv201;
         fphi[i][15] = tuv120;
         fphi[i][16] = tuv021;
         fphi[i][17] = tuv102;
         fphi[i][18] = tuv012;
         fphi[i][19] = tuv111;
      }


      if CONSTEXPR (eq<T, UIND>()) {
         real tuv100_1 = 0;
         real tuv010_1 = 0;
         real tuv001_1 = 0;
         real tuv200_1 = 0;
         real tuv020_1 = 0;
         real tuv002_1 = 0;
         real tuv110_1 = 0;
         real tuv101_1 = 0;
         real tuv011_1 = 0;
         real tuv100_2 = 0;
         real tuv010_2 = 0;
         real tuv001_2 = 0;
         real tuv200_2 = 0;
         real tuv020_2 = 0;
         real tuv002_2 = 0;
         real tuv110_2 = 0;
         real tuv101_2 = 0;
         real tuv011_2 = 0;
         real tuv000 = 0;
         real tuv001 = 0;
         real tuv010 = 0;
         real tuv100 = 0;
         real tuv200 = 0;
         real tuv020 = 0;
         real tuv002 = 0;
         real tuv110 = 0;
         real tuv101 = 0;
         real tuv011 = 0;
         real tuv300 = 0;
         real tuv030 = 0;
         real tuv003 = 0;
         real tuv210 = 0;
         real tuv201 = 0;
         real tuv120 = 0;
         real tuv021 = 0;
         real tuv102 = 0;
         real tuv012 = 0;
         real tuv111 = 0;
         for (int iz = 0; iz < bsorder; ++iz) {
            int zbase = igrid3 + iz;
            zbase -= (zbase >= nfft3 ? nfft3 : 0);
            zbase *= (nfft1 * nfft2);
            real v0 = thetai3[4 * iz];
            real v1 = thetai3[4 * iz + 1];
            real v2 = thetai3[4 * iz + 2];
            real v3 = thetai3[4 * iz + 3];
            real tu00_1 = 0;
            real tu01_1 = 0;
            real tu10_1 = 0;
            real tu20_1 = 0;
            real tu11_1 = 0;
            real tu02_1 = 0;
            real tu00_2 = 0;
            real tu01_2 = 0;
            real tu10_2 = 0;
            real tu20_2 = 0;
            real tu11_2 = 0;
            real tu02_2 = 0;
            real tu00 = 0;
            real tu10 = 0;
            real tu01 = 0;
            real tu20 = 0;
            real tu11 = 0;
            real tu02 = 0;
            real tu30 = 0;
            real tu21 = 0;
            real tu12 = 0;
            real tu03 = 0;
            for (int iy = 0; iy < bsorder; ++iy) {
               int ybase = igrid2 + iy;
               ybase -= (ybase >= nfft2 ? nfft2 : 0);
               ybase *= nfft1;
               real u0 = thetai2[4 * iy];
               real u1 = thetai2[4 * iy + 1];
               real u2 = thetai2[4 * iy + 2];
               real u3 = thetai2[4 * iy + 3];
               real t0_1 = 0;
               real t1_1 = 0;
               real t2_1 = 0;
               real t0_2 = 0;
               real t1_2 = 0;
               real t2_2 = 0;
               real t3 = 0;
               for (int ix = 0; ix < bsorder; ++ix) {
                  int xbase = igrid1 + ix;
                  xbase -= (xbase >= nfft1 ? nfft1 : 0);
                  real tq_1 = qgrid[2 * (xbase + ybase + zbase)];
                  real tq_2 = qgrid[2 * (xbase + ybase + zbase) + 1];
                  t0_1 += tq_1 * thetai1[4 * ix];
                  t1_1 += tq_1 * thetai1[4 * ix + 1];
                  t2_1 += tq_1 * thetai1[4 * ix + 2];
                  t0_2 += tq_2 * thetai1[4 * ix];
                  t1_2 += tq_2 * thetai1[4 * ix + 1];
                  t2_2 += tq_2 * thetai1[4 * ix + 2];
                  t3 += (tq_1 + tq_2) * thetai1[4 * ix + 3];
               }
               tu00_1 += t0_1 * u0;
               tu10_1 += t1_1 * u0;
               tu01_1 += t0_1 * u1;
               tu20_1 += t2_1 * u0;
               tu11_1 += t1_1 * u1;
               tu02_1 += t0_1 * u2;
               tu00_2 += t0_2 * u0;
               tu10_2 += t1_2 * u0;
               tu01_2 += t0_2 * u1;
               tu20_2 += t2_2 * u0;
               tu11_2 += t1_2 * u1;
               tu02_2 += t0_2 * u2;
               real t0 = t0_1 + t0_2;
               real t1 = t1_1 + t1_2;
               real t2 = t2_1 + t2_2;
               tu00 += t0 * u0;
               tu10 += t1 * u0;
               tu01 += t0 * u1;
               tu20 += t2 * u0;
               tu11 += t1 * u1;
               tu02 += t0 * u2;
               tu30 += t3 * u0;
               tu21 += t2 * u1;
               tu12 += t1 * u2;
               tu03 += t0 * u3;
            }
            tuv100_1 += tu10_1 * v0;
            tuv010_1 += tu01_1 * v0;
            tuv001_1 += tu00_1 * v1;
            tuv200_1 += tu20_1 * v0;
            tuv020_1 += tu02_1 * v0;
            tuv002_1 += tu00_1 * v2;
            tuv110_1 += tu11_1 * v0;
            tuv101_1 += tu10_1 * v1;
            tuv011_1 += tu01_1 * v1;
            tuv100_2 += tu10_2 * v0;
            tuv010_2 += tu01_2 * v0;
            tuv001_2 += tu00_2 * v1;
            tuv200_2 += tu20_2 * v0;
            tuv020_2 += tu02_2 * v0;
            tuv002_2 += tu00_2 * v2;
            tuv110_2 += tu11_2 * v0;
            tuv101_2 += tu10_2 * v1;
            tuv011_2 += tu01_2 * v1;
            tuv000 += tu00 * v0;
            tuv100 += tu10 * v0;
            tuv010 += tu01 * v0;
            tuv001 += tu00 * v1;
            tuv200 += tu20 * v0;
            tuv020 += tu02 * v0;
            tuv002 += tu00 * v2;
            tuv110 += tu11 * v0;
            tuv101 += tu10 * v1;
            tuv011 += tu01 * v1;
            tuv300 += tu30 * v0;
            tuv030 += tu03 * v0;
            tuv003 += tu00 * v3;
            tuv210 += tu21 * v0;
            tuv201 += tu20 * v1;
            tuv120 += tu12 * v0;
            tuv021 += tu02 * v1;
            tuv102 += tu10 * v2;
            tuv012 += tu01 * v2;
            tuv111 += tu11 * v1;
         } // end for (iz)
         fdip_phi1[i][0] = 0;
         fdip_phi1[i][1] = tuv100_1;
         fdip_phi1[i][2] = tuv010_1;
         fdip_phi1[i][3] = tuv001_1;
         fdip_phi1[i][4] = tuv200_1;
         fdip_phi1[i][5] = tuv020_1;
         fdip_phi1[i][6] = tuv002_1;
         fdip_phi1[i][7] = tuv110_1;
         fdip_phi1[i][8] = tuv101_1;
         fdip_phi1[i][9] = tuv011_1;
         fdip_phi2[i][0] = 0;
         fdip_phi2[i][1] = tuv100_2;
         fdip_phi2[i][2] = tuv010_2;
         fdip_phi2[i][3] = tuv001_2;
         fdip_phi2[i][4] = tuv200_2;
         fdip_phi2[i][5] = tuv020_2;
         fdip_phi2[i][6] = tuv002_2;
         fdip_phi2[i][7] = tuv110_2;
         fdip_phi2[i][8] = tuv101_2;
         fdip_phi2[i][9] = tuv011_2;
         fdip_sum_phi[i][0] = tuv000;
         fdip_sum_phi[i][1] = tuv100;
         fdip_sum_phi[i][2] = tuv010;
         fdip_sum_phi[i][3] = tuv001;
         fdip_sum_phi[i][4] = tuv200;
         fdip_sum_phi[i][5] = tuv020;
         fdip_sum_phi[i][6] = tuv002;
         fdip_sum_phi[i][7] = tuv110;
         fdip_sum_phi[i][8] = tuv101;
         fdip_sum_phi[i][9] = tuv011;
         fdip_sum_phi[i][10] = tuv300;
         fdip_sum_phi[i][11] = tuv030;
         fdip_sum_phi[i][12] = tuv003;
         fdip_sum_phi[i][13] = tuv210;
         fdip_sum_phi[i][14] = tuv201;
         fdip_sum_phi[i][15] = tuv120;
         fdip_sum_phi[i][16] = tuv021;
         fdip_sum_phi[i][17] = tuv102;
         fdip_sum_phi[i][18] = tuv012;
         fdip_sum_phi[i][19] = tuv111;
      }


      if CONSTEXPR (eq<T, UIND2>()) {
         real tuv100_1 = 0;
         real tuv010_1 = 0;
         real tuv001_1 = 0;
         real tuv200_1 = 0;
         real tuv020_1 = 0;
         real tuv002_1 = 0;
         real tuv110_1 = 0;
         real tuv101_1 = 0;
         real tuv011_1 = 0;
         real tuv100_2 = 0;
         real tuv010_2 = 0;
         real tuv001_2 = 0;
         real tuv200_2 = 0;
         real tuv020_2 = 0;
         real tuv002_2 = 0;
         real tuv110_2 = 0;
         real tuv101_2 = 0;
         real tuv011_2 = 0;
         for (int iz = 0; iz < bsorder; ++iz) {
            int zbase = igrid3 + iz;
            zbase -= (zbase >= nfft3 ? nfft3 : 0);
            zbase *= (nfft1 * nfft2);
            real v0 = thetai3[4 * iz];
            real v1 = thetai3[4 * iz + 1];
            real v2 = thetai3[4 * iz + 2];
            real tu00_1 = 0;
            real tu01_1 = 0;
            real tu10_1 = 0;
            real tu20_1 = 0;
            real tu11_1 = 0;
            real tu02_1 = 0;
            real tu00_2 = 0;
            real tu01_2 = 0;
            real tu10_2 = 0;
            real tu20_2 = 0;
            real tu11_2 = 0;
            real tu02_2 = 0;
            for (int iy = 0; iy < bsorder; ++iy) {
               int ybase = igrid2 + iy;
               ybase -= (ybase >= nfft2 ? nfft2 : 0);
               ybase *= nfft1;
               real u0 = thetai2[4 * iy];
               real u1 = thetai2[4 * iy + 1];
               real u2 = thetai2[4 * iy + 2];
               real t0_1 = 0;
               real t1_1 = 0;
               real t2_1 = 0;
               real t0_2 = 0;
               real t1_2 = 0;
               real t2_2 = 0;
               for (int ix = 0; ix < bsorder; ++ix) {
                  int xbase = igrid1 + ix;
                  xbase -= (xbase >= nfft1 ? nfft1 : 0);
                  real tq_1 = qgrid[2 * (xbase + ybase + zbase)];
                  real tq_2 = qgrid[2 * (xbase + ybase + zbase) + 1];
                  t0_1 += tq_1 * thetai1[4 * ix];
                  t1_1 += tq_1 * thetai1[4 * ix + 1];
                  t2_1 += tq_1 * thetai1[4 * ix + 2];
                  t0_2 += tq_2 * thetai1[4 * ix];
                  t1_2 += tq_2 * thetai1[4 * ix + 1];
                  t2_2 += tq_2 * thetai1[4 * ix + 2];
               }
               tu00_1 += t0_1 * u0;
               tu10_1 += t1_1 * u0;
               tu01_1 += t0_1 * u1;
               tu20_1 += t2_1 * u0;
               tu11_1 += t1_1 * u1;
               tu02_1 += t0_1 * u2;
               tu00_2 += t0_2 * u0;
               tu10_2 += t1_2 * u0;
               tu01_2 += t0_2 * u1;
               tu20_2 += t2_2 * u0;
               tu11_2 += t1_2 * u1;
               tu02_2 += t0_2 * u2;
            }
            tuv100_1 += tu10_1 * v0;
            tuv010_1 += tu01_1 * v0;
            tuv001_1 += tu00_1 * v1;
            tuv200_1 += tu20_1 * v0;
            tuv020_1 += tu02_1 * v0;
            tuv002_1 += tu00_1 * v2;
            tuv110_1 += tu11_1 * v0;
            tuv101_1 += tu10_1 * v1;
            tuv011_1 += tu01_1 * v1;
            tuv100_2 += tu10_2 * v0;
            tuv010_2 += tu01_2 * v0;
            tuv001_2 += tu00_2 * v1;
            tuv200_2 += tu20_2 * v0;
            tuv020_2 += tu02_2 * v0;
            tuv002_2 += tu00_2 * v2;
            tuv110_2 += tu11_2 * v0;
            tuv101_2 += tu10_2 * v1;
            tuv011_2 += tu01_2 * v1;
         } // end for (iz)
         fdip_phi1[i][0] = 0;
         fdip_phi1[i][1] = tuv100_1;
         fdip_phi1[i][2] = tuv010_1;
         fdip_phi1[i][3] = tuv001_1;
         fdip_phi1[i][4] = tuv200_1;
         fdip_phi1[i][5] = tuv020_1;
         fdip_phi1[i][6] = tuv002_1;
         fdip_phi1[i][7] = tuv110_1;
         fdip_phi1[i][8] = tuv101_1;
         fdip_phi1[i][9] = tuv011_1;
         fdip_phi2[i][0] = 0;
         fdip_phi2[i][1] = tuv100_2;
         fdip_phi2[i][2] = tuv010_2;
         fdip_phi2[i][3] = tuv001_2;
         fdip_phi2[i][4] = tuv200_2;
         fdip_phi2[i][5] = tuv020_2;
         fdip_phi2[i][6] = tuv002_2;
         fdip_phi2[i][7] = tuv110_2;
         fdip_phi2[i][8] = tuv101_2;
         fdip_phi2[i][9] = tuv011_2;
      }
   }
}


void fphi_mpole_cu(PMEUnit pme_u, real (*fphi)[20])
{
   auto& st = *pme_u;
   int n1 = st.nfft1;
   int n2 = st.nfft2;
   int n3 = st.nfft3;


   auto ker = fphi_get_cu<MPOLE, 5>;
   launch_k2s(nonblk, PME_BLOCKDIM, n, ker, n, n1, n2, n3, x, y, z, (real*)fphi,
              nullptr, nullptr, st.qgrid, recipa, recipb, recipc);
}


void fphi_uind_cu(PMEUnit pme_u, real (*fdip_phi1)[10], real (*fdip_phi2)[10],
                  real (*fdip_sum_phi)[20])
{
   auto& st = *pme_u;
   int n1 = st.nfft1;
   int n2 = st.nfft2;
   int n3 = st.nfft3;


   auto ker = fphi_get_cu<UIND, 5>;
   launch_k2s(nonblk, PME_BLOCKDIM, n, ker, n, n1, n2, n3, x, y, z,
              (real*)fdip_phi1, (real*)fdip_phi2, (real*)fdip_sum_phi, st.qgrid,
              recipa, recipb, recipc);
}


void fphi_uind2_cu(PMEUnit pme_u, real (*fdip_phi1)[10], real (*fdip_phi2)[10])
{
   auto& st = *pme_u;
   int n1 = st.nfft1;
   int n2 = st.nfft2;
   int n3 = st.nfft3;


   auto ker = fphi_get_cu<UIND2, 5>;
   launch_k2s(nonblk, PME_BLOCKDIM, n, ker, n, n1, n2, n3, x, y, z,
              (real*)fdip_phi1, (real*)fdip_phi2, nullptr, st.qgrid, recipa,
              recipb, recipc);
}


void pme_cuda_func_config()
{
   // grid

   auto grid_mpolek = grid_put_cu1<MPOLE, 5>;
   check_rt(cudaFuncSetCacheConfig(grid_mpolek, cudaFuncCachePreferNone));

   auto grid_uindk = grid_put_cu1<UIND, 5>;
   check_rt(cudaFuncSetCacheConfig(grid_uindk, cudaFuncCachePreferNone));

   // fphi

   auto fphi_mpole = fphi_get_cu<MPOLE, 5>;
   check_rt(cudaFuncSetCacheConfig(fphi_mpole, cudaFuncCachePreferL1));

   auto fphi_uind = fphi_get_cu<UIND, 5>;
   check_rt(cudaFuncSetCacheConfig(fphi_uind, cudaFuncCachePreferL1));

   auto fphi_uind2 = fphi_get_cu<UIND2, 5>;
   check_rt(cudaFuncSetCacheConfig(fphi_uind2, cudaFuncCachePreferL1));
}
TINKER_NAMESPACE_END
