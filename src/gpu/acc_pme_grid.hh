#ifdef TINKER_SRC_GPU_ACC_PMESTUFF_IMPL_

TINKER_NAMESPACE_BEGIN
namespace gpu {
// see also subroutine bsplgen in pmestuf.f
#pragma acc routine seq
template <int LEVEL>
static inline void bsplgen(real w, real* __restrict__ thetai,
                           real* __restrict__ _bsbuild, int bsorder) {
  static_assert(LEVEL == 2 || LEVEL == 4, "");

  real_allocatable bsbuild(_bsbuild, bsorder, bsorder);

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

  // get coefficients for the B-spline first derivative

  int k = bsorder - 1;
  bsbuild(k, bsorder) = bsbuild(k, bsorder - 1);
  for (int i = bsorder - 1; i >= 2; --i) {
    bsbuild(k, i) = bsbuild(k, i - 1) - bsbuild(k, i);
  }
  bsbuild(k, 1) = -bsbuild(k, 1);

  if_constexpr(LEVEL == 4) {

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
    // LEVEL = 2; j = 1, 2
    thetai[4 * (i - 1)] = bsbuild(bsorder, i);
    thetai[4 * (i - 1) + 1] = bsbuild(bsorder - 1, i);
    if_constexpr(LEVEL == 4) {
      // LEVEL = 4; j = 3, 4
      thetai[4 * (i - 1) + 2] = bsbuild(bsorder - 2, i);
      thetai[4 * (i - 1) + 3] = bsbuild(bsorder - 3, i);
    }
  }
}

void grid_mpole(real (*fmp)[10]) {
  pme_st& st = pme_obj(epme_unit);
  pme_st* dptr = pme_deviceptr(epme_unit);

  const int nfft1 = st.nfft1;
  const int nfft2 = st.nfft2;
  const int nfft3 = st.nfft3;
  const int bsorder = st.bsorder;

  const int order4 = 4 * bsorder;
  static std::vector<real> t1buf, t2buf, t3buf;
  t1buf.resize(order4);
  t2buf.resize(order4);
  t3buf.resize(order4);
  real *thetai1 = t1buf.data(), *thetai2 = t2buf.data(),
       *thetai3 = t3buf.data();

  const int bso2 = bsorder * bsorder;
  static std::vector<real> bsbuildbuf;
  bsbuildbuf.resize(bso2);
  real* bsbuild = bsbuildbuf.data();

  #pragma acc parallel loop independent deviceptr(dptr)
  for (int i = 0; i < 2 * nfft1 * nfft2 * nfft3; ++i) {
    dptr->qgrid[i] = 0;
  } // for (int i)

  #pragma acc parallel loop independent deviceptr(fmp,x,y,z,box,dptr)\
              private(bsbuild[0:bso2],\
              thetai1[0:order4],thetai2[0:order4],thetai3[0:order4])
  for (int i = 0; i < n; ++i) {
    fmat_real3 recip(box->recip);
    real xi = x[i];
    real yi = y[i];
    real zi = z[i];

    // map fractional coordinate w from [-0.5 + k, 0.5 + k) to [0, 1)
    // w -> (w + 0.5) - FLOOR(w + 0.5)
    // see also subroutine bspline_fill in pmestuf.f

    real w1 = xi * recip(1, 1) + yi * recip(2, 1) + zi * recip(3, 1);
    w1 = w1 + 0.5f - REAL_FLOOR(w1 + 0.5f);
    real fr1 = nfft1 * w1;
    int igrid1 = REAL_FLOOR(fr1);
    w1 = fr1 - igrid1;
    bsplgen<4>(w1, thetai1, bsbuild, bsorder);

    real w2 = xi * recip(1, 2) + yi * recip(2, 2) + zi * recip(3, 2);
    w2 = w2 + 0.5f - REAL_FLOOR(w2 + 0.5f);
    real fr2 = nfft2 * w2;
    int igrid2 = REAL_FLOOR(fr2);
    w2 = fr2 - igrid2;
    bsplgen<4>(w2, thetai2, bsbuild, bsorder);

    real w3 = xi * recip(1, 3) + yi * recip(2, 3) + zi * recip(3, 3);
    w3 = w3 + 0.5f - REAL_FLOOR(w3 + 0.5f);
    real fr3 = nfft3 * w3;
    int igrid3 = REAL_FLOOR(fr3);
    w3 = fr3 - igrid3;
    bsplgen<4>(w3, thetai3, bsbuild, bsorder);

    igrid1 = igrid1 - bsorder + 1;
    igrid2 = igrid2 - bsorder + 1;
    igrid3 = igrid3 - bsorder + 1;
    igrid1 += (igrid1 < 0 ? nfft1 : 0);
    igrid2 += (igrid2 < 0 ? nfft2 : 0);
    igrid3 += (igrid3 < 0 ? nfft3 : 0);

    #pragma acc loop independent
    for (int iz = 0; iz < bsorder; ++iz) {
      int zbase = igrid3 + iz;
      zbase -= (zbase >= nfft3 ? nfft3 : 0);
      zbase *= (nfft1 * nfft2);
      real v0 = thetai3[4 * iz];
      real v1 = thetai3[4 * iz + 1];
      real v2 = thetai3[4 * iz + 2];
      #pragma acc loop independent
      for (int iy = 0; iy < bsorder; ++iy) {
        int ybase = igrid2 + iy;
        ybase -= (ybase >= nfft2 ? nfft2 : 0);
        ybase *= nfft1;
        real u0 = thetai2[4 * iy];
        real u1 = thetai2[4 * iy + 1];
        real u2 = thetai2[4 * iy + 2];
        // fmp: 0, x, y, z, xx, yy, zz, xy, xz, yz
        //      1, 2, 3, 4,  5,  6,  7,  8,  9, 10
        real term0 = fmp[i][mpl_pme_0] * u0 * v0 + fmp[i][mpl_pme_y] * u1 * v0 +
            fmp[i][mpl_pme_z] * u0 * v1 + fmp[i][mpl_pme_yy] * u2 * v0 +
            fmp[i][mpl_pme_zz] * u0 * v2 + fmp[i][mpl_pme_yz] * u1 * v1;
        real term1 = fmp[i][mpl_pme_x] * u0 * v0 +
            fmp[i][mpl_pme_xy] * u1 * v0 + fmp[i][mpl_pme_xz] * u0 * v1;
        real term2 = fmp[i][mpl_pme_xx] * u0 * v0;
        #pragma acc loop independent
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
    }
  } // for (int i)
}

void fphi_mpole(real (*fphi)[20]) {
  pme_st& st = pme_obj(epme_unit);
  pme_st* dptr = pme_deviceptr(epme_unit);

  const int nfft1 = st.nfft1;
  const int nfft2 = st.nfft2;
  const int nfft3 = st.nfft3;
  const int bsorder = st.bsorder;

  const int order4 = 4 * bsorder;
  static std::vector<real> t1buf, t2buf, t3buf;
  t1buf.resize(order4);
  t2buf.resize(order4);
  t3buf.resize(order4);
  real *thetai1 = t1buf.data(), *thetai2 = t2buf.data(),
       *thetai3 = t3buf.data();

  const int bso2 = bsorder * bsorder;
  static std::vector<real> bsbuildbuf;
  bsbuildbuf.resize(bso2);
  real* bsbuild = bsbuildbuf.data();

  #pragma acc parallel loop independent deviceptr(fphi,x,y,z,box,dptr)\
              private(bsbuild[0:bso2],\
              thetai1[0:order4],thetai2[0:order4],thetai3[0:order4])
  for (int i = 0; i < n; ++i) {
    fmat_real3 recip(box->recip);
    real xi = x[i];
    real yi = y[i];
    real zi = z[i];

    // map fractional coordinate w from [-0.5 + k, 0.5 + k) to [0, 1)
    // w -> (w + 0.5) - FLOOR(w + 0.5)
    // see also subroutine bspline_fill in pmestuf.f

    real w1 = xi * recip(1, 1) + yi * recip(2, 1) + zi * recip(3, 1);
    w1 = w1 + 0.5f - REAL_FLOOR(w1 + 0.5f);
    real fr1 = nfft1 * w1;
    int igrid1 = REAL_FLOOR(fr1);
    w1 = fr1 - igrid1;
    bsplgen<4>(w1, thetai1, bsbuild, bsorder);

    real w2 = xi * recip(1, 2) + yi * recip(2, 2) + zi * recip(3, 2);
    w2 = w2 + 0.5f - REAL_FLOOR(w2 + 0.5f);
    real fr2 = nfft2 * w2;
    int igrid2 = REAL_FLOOR(fr2);
    w2 = fr2 - igrid2;
    bsplgen<4>(w2, thetai2, bsbuild, bsorder);

    real w3 = xi * recip(1, 3) + yi * recip(2, 3) + zi * recip(3, 3);
    w3 = w3 + 0.5f - REAL_FLOOR(w3 + 0.5f);
    real fr3 = nfft3 * w3;
    int igrid3 = REAL_FLOOR(fr3);
    w3 = fr3 - igrid3;
    bsplgen<4>(w3, thetai3, bsbuild, bsorder);

    igrid1 = igrid1 - bsorder + 1;
    igrid2 = igrid2 - bsorder + 1;
    igrid3 = igrid3 - bsorder + 1;
    igrid1 += (igrid1 < 0 ? nfft1 : 0);
    igrid2 += (igrid2 < 0 ? nfft2 : 0);
    igrid3 += (igrid3 < 0 ? nfft3 : 0);

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
    #pragma acc loop independent reduction(+:\
                tuv000,tuv001,tuv010,tuv100,tuv200,\
                tuv020,tuv002,tuv110,tuv101,tuv011,\
                tuv300,tuv030,tuv003,tuv210,tuv201,\
                tuv120,tuv021,tuv102,tuv012,tuv111)
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
      #pragma acc loop independent reduction(+:\
                  tu00,tu10,tu01,tu20,tu11,tu02,tu30,tu21,tu12,tu03)
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
        #pragma acc loop independent reduction(+:t0,t1,t2,t3)
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
  }
}
}
TINKER_NAMESPACE_END

#endif
