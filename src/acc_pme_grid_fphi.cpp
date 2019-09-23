#include "acc_mathfunc.h"
#include "box.h"
#include "e_mpole.h"
#include "md.h"
#include "pme.h"

TINKER_NAMESPACE_BEGIN
static constexpr int MAX_BSORDER = 5;

// see also subroutine bsplgen in pmestuf.f
#pragma acc routine seq
template <int LEVEL>
static inline void bsplgen(real w, real* __restrict__ thetai, int bsorder) {

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

#define bsbuild(j, i) bsbuild_[((i)-1) * bsorder + (j)-1]

  real bsbuild_[MAX_BSORDER * MAX_BSORDER];

  // initialization to get to 2nd order recursion

  bsbuild(2, 2) = w;
  bsbuild(2, 1) = 1 - w;

  // perform one pass to get to 3rd order recursion

  bsbuild(3, 3) = 0.5f * w * bsbuild(2, 2);
  bsbuild(3, 2) = 0.5f * ((1 + w) * bsbuild(2, 1) + (2 - w) * bsbuild(2, 2));
  bsbuild(3, 1) = 0.5f * (1 - w) * bsbuild(2, 1);

  // compute standard B-spline recursion to desired order

  #pragma acc loop seq
  for (int i = 4; i <= bsorder; ++i) {
    int k = i - 1;
    real denom = REAL_RECIP(k);
    bsbuild(i, i) = denom * w * bsbuild(k, k);
    #pragma acc loop seq
    for (int j = 1; j <= i - 2; j++) {
      bsbuild(i, i - j) = denom *
          ((w + j) * bsbuild(k, i - j - 1) + (i - j - w) * bsbuild(k, i - j));
    }
    bsbuild(i, 1) = denom * (1 - w) * bsbuild(k, 1);
  }

  int k;

  if_constexpr(LEVEL >= 2) {

    // get coefficients for the B-spline first derivative

    k = bsorder - 1;
    bsbuild(k, bsorder) = bsbuild(k, bsorder - 1);
    #pragma acc loop seq
    for (int i = bsorder - 1; i >= 2; --i) {
      bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
    }
    bsbuild(k, 1) = -bsbuild(k, 1);
  }

  if_constexpr(LEVEL >= 3) {

    // get coefficients for the B-spline second derivative

    k = bsorder - 2;
    bsbuild(k, bsorder - 1) = bsbuild(k, bsorder - 2);
    #pragma acc loop seq
    for (int i = bsorder - 2; i >= 2; --i) {
      bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
    }
    bsbuild(k, 1) = -bsbuild(k, 1);
    bsbuild(k, bsorder) = bsbuild(k, bsorder - 1);
    #pragma acc loop seq
    for (int i = bsorder - 1; i >= 2; --i) {
      bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
    }
    bsbuild(k, 1) = -bsbuild(k, 1);
  }

  if_constexpr(LEVEL == 4) {

    // get coefficients for the B-spline third derivative

    k = bsorder - 3;
    bsbuild(k, bsorder - 2) = bsbuild(k, bsorder - 3);
    #pragma acc loop seq
    for (int i = bsorder - 3; i >= 2; --i) {
      bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
    }
    bsbuild(k, 1) = -bsbuild(k, 1);
    bsbuild(k, bsorder - 1) = bsbuild(k, bsorder - 2);
    #pragma acc loop seq
    for (int i = bsorder - 2; i >= 2; --i) {
      bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
    }
    bsbuild(k, 1) = -bsbuild(k, 1);
    bsbuild(k, bsorder) = bsbuild(k, bsorder - 1);
    #pragma acc loop seq
    for (int i = bsorder - 1; i >= 2; --i)
      bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
    bsbuild(k, 1) = -bsbuild(k, 1);
  }

  // copy coefficients from temporary to permanent storage

  #pragma acc loop seq
  for (int i = 1; i <= bsorder; ++i) {
    #pragma acc loop seq
    for (int j = 1; j <= LEVEL; ++j) {
      thetai[4 * (i - 1) + (j - 1)] = bsbuild(bsorder - j + 1, i);
    }
  }

#undef bsbuild
}

enum { PCHG_GRID = 1, MPOLE_GRID, UIND_GRID, UIND_GRID_FPHI2, DISP_GRID };

template <int WHAT>
void grid_tmpl(PMEUnit pme_u, real* optional1, real* optional2) {
  auto& st = *pme_u;
  auto* dptr = pme_u.deviceptr();

  MAYBE_UNUSED real* pchg;
  if_constexpr(WHAT == PCHG_GRID || WHAT == DISP_GRID) { pchg = optional1; }

  MAYBE_UNUSED real(*fmp)[10];
  if_constexpr(WHAT == MPOLE_GRID) {
    fmp = reinterpret_cast<real(*)[10]>(optional1);
  }

  MAYBE_UNUSED real(*fuind)[3];
  MAYBE_UNUSED real(*fuinp)[3];
  if_constexpr(WHAT == UIND_GRID) {
    fuind = reinterpret_cast<real(*)[3]>(optional1);
    fuinp = reinterpret_cast<real(*)[3]>(optional2);
  }

  const int nfft1 = st.nfft1;
  const int nfft2 = st.nfft2;
  const int nfft3 = st.nfft3;
  const int bsorder = st.bsorder;
  assert(bsorder <= MAX_BSORDER);

  device_array::zero(2 * nfft1 * nfft2 * nfft3, st.qgrid);

  #pragma acc parallel loop independent\
              deviceptr(pchg,fmp,fuind,fuinp,x,y,z,box,dptr)
  for (int i = 0; i < n; ++i) {
    real thetai1[4 * MAX_BSORDER], thetai2[4 * MAX_BSORDER],
        thetai3[4 * MAX_BSORDER];

    real xi = x[i];
    real yi = y[i];
    real zi = z[i];

    // map fractional coordinate w from [-0.5 + k, 0.5 + k) to [0, 1)
    // w -> (w + 0.5) - FLOOR(w + 0.5)
    // see also subroutine bspline_fill in pmestuf.f

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

    if_constexpr(WHAT == PCHG_GRID || WHAT == DISP_GRID) {
      bsplgen<1>(w1, thetai1, bsorder);
      bsplgen<1>(w2, thetai2, bsorder);
      bsplgen<1>(w3, thetai3, bsorder);
    }

    if_constexpr(WHAT == MPOLE_GRID) {
      bsplgen<3>(w1, thetai1, bsorder);
      bsplgen<3>(w2, thetai2, bsorder);
      bsplgen<3>(w3, thetai3, bsorder);
    }

    if_constexpr(WHAT == UIND_GRID) {
      bsplgen<2>(w1, thetai1, bsorder);
      bsplgen<2>(w2, thetai2, bsorder);
      bsplgen<2>(w3, thetai3, bsorder);
    }

    igrid1 = igrid1 - bsorder + 1;
    igrid2 = igrid2 - bsorder + 1;
    igrid3 = igrid3 - bsorder + 1;
    igrid1 += (igrid1 < 0 ? nfft1 : 0);
    igrid2 += (igrid2 < 0 ? nfft2 : 0);
    igrid3 += (igrid3 < 0 ? nfft3 : 0);

    if_constexpr(WHAT == PCHG_GRID || WHAT == DISP_GRID) {
      #pragma acc loop seq
      for (int iz = 0; iz < bsorder; ++iz) {
        int zbase = igrid3 + iz;
        zbase -= (zbase >= nfft3 ? nfft3 : 0);
        zbase *= (nfft1 * nfft2);
        real v0 = thetai3[4 * iz] * pchg[i];
        #pragma acc loop seq
        for (int iy = 0; iy < bsorder; ++iy) {
          int ybase = igrid2 + iy;
          ybase -= (ybase >= nfft2 ? nfft2 : 0);
          ybase *= nfft1;
          real u0 = thetai2[4 * iy];
          real term = v0 * u0;
          #pragma acc loop seq
          for (int ix = 0; ix < bsorder; ++ix) {
            int xbase = igrid1 + ix;
            xbase -= (xbase >= nfft1 ? nfft1 : 0);
            int index = xbase + ybase + zbase;
            real t0 = thetai1[4 * ix];
            #pragma acc atomic update
            dptr->qgrid[2 * index] += term * t0;
          }
        }
      } // end for (int iz)
    }   // end if (grid_pchg || grid_disp)

    if_constexpr(WHAT == MPOLE_GRID) {
      #pragma acc loop seq
      for (int iz = 0; iz < bsorder; ++iz) {
        int zbase = igrid3 + iz;
        zbase -= (zbase >= nfft3 ? nfft3 : 0);
        zbase *= (nfft1 * nfft2);
        real v0 = thetai3[4 * iz];
        real v1 = thetai3[4 * iz + 1];
        real v2 = thetai3[4 * iz + 2];
        #pragma acc loop seq
        for (int iy = 0; iy < bsorder; ++iy) {
          int ybase = igrid2 + iy;
          ybase -= (ybase >= nfft2 ? nfft2 : 0);
          ybase *= nfft1;
          real u0 = thetai2[4 * iy];
          real u1 = thetai2[4 * iy + 1];
          real u2 = thetai2[4 * iy + 2];
          // fmp: 0, x, y, z, xx, yy, zz, xy, xz, yz
          //      1, 2, 3, 4,  5,  6,  7,  8,  9, 10
          real term0 = fmp[i][mpl_pme_0] * u0 * v0 +
              fmp[i][mpl_pme_y] * u1 * v0 + fmp[i][mpl_pme_z] * u0 * v1 +
              fmp[i][mpl_pme_yy] * u2 * v0 + fmp[i][mpl_pme_zz] * u0 * v2 +
              fmp[i][mpl_pme_yz] * u1 * v1;
          real term1 = fmp[i][mpl_pme_x] * u0 * v0 +
              fmp[i][mpl_pme_xy] * u1 * v0 + fmp[i][mpl_pme_xz] * u0 * v1;
          real term2 = fmp[i][mpl_pme_xx] * u0 * v0;
          #pragma acc loop seq
          for (int ix = 0; ix < bsorder; ++ix) {
            int xbase = igrid1 + ix;
            xbase -= (xbase >= nfft1 ? nfft1 : 0);
            int index = xbase + ybase + zbase;
            real t0 = thetai1[4 * ix];
            real t1 = thetai1[4 * ix + 1];
            real t2 = thetai1[4 * ix + 2];
            #pragma acc atomic update
            dptr->qgrid[2 * index] += term0 * t0 + term1 * t1 + term2 * t2;
          }
        }
      } // end for (int iz)
    }   // end if (grid_mpole)

    if_constexpr(WHAT == UIND_GRID) {
      #pragma acc loop seq
      for (int iz = 0; iz < bsorder; ++iz) {
        int zbase = igrid3 + iz;
        zbase -= (zbase >= nfft3 ? nfft3 : 0);
        zbase *= (nfft1 * nfft2);
        real v0 = thetai3[4 * iz];
        real v1 = thetai3[4 * iz + 1];
        #pragma acc loop seq
        for (int iy = 0; iy < bsorder; ++iy) {
          int ybase = igrid2 + iy;
          ybase -= (ybase >= nfft2 ? nfft2 : 0);
          ybase *= nfft1;
          real u0 = thetai2[4 * iy];
          real u1 = thetai2[4 * iy + 1];
          real term01 = fuind[i][1] * u1 * v0 + fuind[i][2] * u0 * v1;
          real term11 = fuind[i][0] * u0 * v0;
          real term02 = fuinp[i][1] * u1 * v0 + fuinp[i][2] * u0 * v1;
          real term12 = fuinp[i][0] * u0 * v0;
          #pragma acc loop seq
          for (int ix = 0; ix < bsorder; ++ix) {
            int xbase = igrid1 + ix;
            xbase -= (xbase >= nfft1 ? nfft1 : 0);
            int index = xbase + ybase + zbase;
            real t0 = thetai1[4 * ix];
            real t1 = thetai1[4 * ix + 1];
            #pragma acc atomic update
            dptr->qgrid[2 * index] += (term01 * t0 + term11 * t1);
            #pragma acc atomic update
            dptr->qgrid[2 * index + 1] += (term02 * t0 + term12 * t1);
          }
        }
      } // end for (int iz)
    }   // end if (grid_uind)
  }     // for (int i)
}

void grid_mpole(PMEUnit pme_u, real (*fmp)[10]) {
  real* opt1 = reinterpret_cast<real*>(fmp);
  grid_tmpl<MPOLE_GRID>(pme_u, opt1, nullptr);
}

void grid_uind(PMEUnit pme_u, real (*fuind)[3], real (*fuinp)[3]) {
  real* opt1 = reinterpret_cast<real*>(fuind);
  real* opt2 = reinterpret_cast<real*>(fuinp);
  grid_tmpl<UIND_GRID>(pme_u, opt1, opt2);
}

template <int WHAT>
void fphi_tmpl(PMEUnit pme_u, real* opt1, real* opt2, real* opt3) {
  auto& st = *pme_u;
  auto* dptr = pme_u.deviceptr();

  MAYBE_UNUSED real(*fphi)[20];
  if_constexpr(WHAT == MPOLE_GRID) {
    fphi = reinterpret_cast<real(*)[20]>(opt1);
  }

  MAYBE_UNUSED real(*fdip_phi1)[10];
  MAYBE_UNUSED real(*fdip_phi2)[10];
  MAYBE_UNUSED real(*fdip_sum_phi)[20];
  if_constexpr(WHAT == UIND_GRID) {
    fdip_phi1 = reinterpret_cast<real(*)[10]>(opt1);
    fdip_phi2 = reinterpret_cast<real(*)[10]>(opt2);
    fdip_sum_phi = reinterpret_cast<real(*)[20]>(opt3);
  }
  if_constexpr(WHAT == UIND_GRID_FPHI2) {
    fdip_phi1 = reinterpret_cast<real(*)[10]>(opt1);
    fdip_phi2 = reinterpret_cast<real(*)[10]>(opt2);
  }

  const int nfft1 = st.nfft1;
  const int nfft2 = st.nfft2;
  const int nfft3 = st.nfft3;
  const int bsorder = st.bsorder;
  assert(bsorder <= MAX_BSORDER);

  #pragma acc parallel loop independent\
              deviceptr(fphi,fdip_phi1,fdip_phi2,fdip_sum_phi,\
              x,y,z,box,dptr)
  for (int i = 0; i < n; ++i) {
    real thetai1[4 * MAX_BSORDER], thetai2[4 * MAX_BSORDER],
        thetai3[4 * MAX_BSORDER];

    real xi = x[i];
    real yi = y[i];
    real zi = z[i];

    // map fractional coordinate w from [-0.5 + k, 0.5 + k) to [0, 1)
    // w -> (w + 0.5) - FLOOR(w + 0.5)
    // see also subroutine bspline_fill in pmestuf.f

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

    if_constexpr(WHAT == MPOLE_GRID || WHAT == UIND_GRID ||
                 WHAT == UIND_GRID_FPHI2) {
      bsplgen<4>(w1, thetai1, bsorder);
      bsplgen<4>(w2, thetai2, bsorder);
      bsplgen<4>(w3, thetai3, bsorder);
    }

    igrid1 = igrid1 - bsorder + 1;
    igrid2 = igrid2 - bsorder + 1;
    igrid3 = igrid3 - bsorder + 1;
    igrid1 += (igrid1 < 0 ? nfft1 : 0);
    igrid2 += (igrid2 < 0 ? nfft2 : 0);
    igrid3 += (igrid3 < 0 ? nfft3 : 0);

    if_constexpr(WHAT == MPOLE_GRID) {
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
      #pragma acc loop seq
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
        #pragma acc loop seq
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
          #pragma acc loop seq
          for (int ix = 0; ix < bsorder; ++ix) {
            int xbase = igrid1 + ix;
            xbase -= (xbase >= nfft1 ? nfft1 : 0);
            real tq = dptr->qgrid[2 * (xbase + ybase + zbase)];
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
    } // end if (fphi_mpole)

    if_constexpr(WHAT == UIND_GRID || WHAT == UIND_GRID_FPHI2) {
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
      #pragma acc loop seq
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
        #pragma acc loop seq
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
          #pragma acc loop seq
          for (int ix = 0; ix < bsorder; ++ix) {
            int xbase = igrid1 + ix;
            xbase -= (xbase >= nfft1 ? nfft1 : 0);
            real tq_1 = dptr->qgrid[2 * (xbase + ybase + zbase)];
            real tq_2 = dptr->qgrid[2 * (xbase + ybase + zbase) + 1];
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

      if_constexpr(WHAT == UIND_GRID) {
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
    } // end if (fphi_uind)
  }   // end for (int i)
}

void fphi_mpole(PMEUnit pme_u, real (*fphi)[20]) {
  real* opt1 = reinterpret_cast<real*>(fphi);
  fphi_tmpl<MPOLE_GRID>(pme_u, opt1, nullptr, nullptr);
}

void fphi_uind(PMEUnit pme_u, real (*fdip_phi1)[10], real (*fdip_phi2)[10],
               real (*fdip_sum_phi)[20]) {
  real* opt1 = reinterpret_cast<real*>(fdip_phi1);
  real* opt2 = reinterpret_cast<real*>(fdip_phi2);
  real* opt3 = reinterpret_cast<real*>(fdip_sum_phi);
  fphi_tmpl<UIND_GRID>(pme_u, opt1, opt2, opt3);
}

void fphi_uind2(PMEUnit pme_u, real (*fdip_phi1)[10], real (*fdip_phi2)[10]) {
  real* opt1 = reinterpret_cast<real*>(fdip_phi1);
  real* opt2 = reinterpret_cast<real*>(fdip_phi2);
  fphi_tmpl<UIND_GRID_FPHI2>(pme_u, opt1, opt2, nullptr);
}
TINKER_NAMESPACE_END
