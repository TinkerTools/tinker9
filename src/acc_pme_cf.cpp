#include "box.h"
#include "gpu/e_mpole.h"
#include "md.h"
#include "pme.h"

/**
 * @file
 * Conversion between Cartesian and fractional coordinates.
 */

TINKER_NAMESPACE_BEGIN
void rpole_to_cmp() {

  // copy multipole moments and coordinates to local storage

  #pragma acc parallel loop independent deviceptr(rpole,cmp)
  for (int i = 0; i < n; ++i) {
    cmp[i][0] = rpole[i][mpl_pme_0];
    cmp[i][1] = rpole[i][mpl_pme_x];
    cmp[i][2] = rpole[i][mpl_pme_y];
    cmp[i][3] = rpole[i][mpl_pme_z];
    cmp[i][4] = rpole[i][mpl_pme_xx];
    cmp[i][5] = rpole[i][mpl_pme_yy];
    cmp[i][6] = rpole[i][mpl_pme_zz];
    cmp[i][7] = 2 * rpole[i][mpl_pme_xy];
    cmp[i][8] = 2 * rpole[i][mpl_pme_xz];
    cmp[i][9] = 2 * rpole[i][mpl_pme_yz];
  }
}

void cmp_to_fmp(PMEUnit pme_u, const real (*cmp)[10], real (*fmp)[10]) {
  auto& st = pme_u.obj();
  int nfft1 = st.nfft1;
  int nfft2 = st.nfft2;
  int nfft3 = st.nfft3;

  // data qi1  / 1, 2, 3, 1, 1, 2 /
  // data qi2  / 1, 2, 3, 2, 3, 3 /
  constexpr int qi1[] = {0, 1, 2, 0, 0, 1};
  constexpr int qi2[] = {0, 1, 2, 1, 2, 2};

  real ctf3[3][3];
  real ctf6[6][6];
  real a[3][3];

  #pragma acc parallel loop deviceptr(box,cmp,fmp)\
              private(ctf3[0:3][0:3],ctf6[0:6][0:6],a[0:3][0:3])
  for (int iatom = 0; iatom < n; ++iatom) {

    // see also subroutine cart_to_frac in pmestuf.f
    // set the reciprocal vector transformation matrix

    a[0][0] = nfft1 * box->recip[0][0];
    a[0][1] = nfft2 * box->recip[1][0];
    a[0][2] = nfft3 * box->recip[2][0];
    a[1][0] = nfft1 * box->recip[0][1];
    a[1][1] = nfft2 * box->recip[1][1];
    a[1][2] = nfft3 * box->recip[2][1];
    a[2][0] = nfft1 * box->recip[0][2];
    a[2][1] = nfft2 * box->recip[1][2];
    a[2][2] = nfft3 * box->recip[2][2];

    // get the Cartesian to fractional conversion matrix

    #pragma acc loop independent collapse(2)
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        ctf3[j][i] = a[j][i];
      }
    }
    #pragma acc loop independent
    for (int i1 = 0; i1 < 3; ++i1) {
      int k = qi1[i1];
      #pragma acc loop independent
      for (int i2 = 0; i2 < 6; ++i2) {
        int i = qi1[i2];
        int j = qi2[i2];
        ctf6[i2][i1] = a[i][k] * a[j][k];
      }
    }
    #pragma acc loop independent
    for (int i1 = 3; i1 < 6; ++i1) {
      int k = qi1[i1];
      int m = qi2[i1];
      #pragma acc loop independent
      for (int i2 = 0; i2 < 6; ++i2) {
        int i = qi1[i2];
        int j = qi2[i2];
        ctf6[i2][i1] = a[i][k] * a[j][m] + a[j][k] * a[i][m];
      }
    }

    // apply the transformation to get the fractional multipoles

    fmp[iatom][0] = cmp[iatom][0];
    #pragma acc loop independent
    for (int j = 1; j < 10; ++j) {
      fmp[iatom][j] = 0;
    }
    #pragma acc loop seq collapse(2)
    for (int j = 1; j < 4; ++j) {
      for (int k = 1; k < 4; ++k) {
        fmp[iatom][j] += ctf3[k - 1][j - 1] * cmp[iatom][k];
      }
    }
    #pragma acc loop seq collapse(2)
    for (int j = 4; j < 10; ++j) {
      for (int k = 4; k < 10; ++k) {
        fmp[iatom][j] += ctf6[k - 4][j - 4] * cmp[iatom][k];
      }
    }
  }
}

void cuind_to_fuind(PMEUnit pme_u, const real (*cind)[3], const real (*cinp)[3],
                    real (*fuind)[3], real (*fuinp)[3]) {
  auto& st = pme_u.obj();
  int nfft1 = st.nfft1;
  int nfft2 = st.nfft2;
  int nfft3 = st.nfft3;

  real a[3][3];

  #pragma acc parallel loop independent deviceptr(box,\
              cind,cinp,fuind,fuinp)\
              private(a[0:3][0:3])
  for (int i = 0; i < n; ++i) {
    a[0][0] = nfft1 * box->recip[0][0];
    a[0][1] = nfft2 * box->recip[1][0];
    a[0][2] = nfft3 * box->recip[2][0];
    a[1][0] = nfft1 * box->recip[0][1];
    a[1][1] = nfft2 * box->recip[1][1];
    a[1][2] = nfft3 * box->recip[2][1];
    a[2][0] = nfft1 * box->recip[0][2];
    a[2][1] = nfft2 * box->recip[1][2];
    a[2][2] = nfft3 * box->recip[2][2];

    #pragma acc loop independent
    for (int j = 0; j < 3; ++j) {
      fuind[i][j] =
          a[0][j] * cind[i][0] + a[1][j] * cind[i][1] + a[2][j] * cind[i][2];
      fuinp[i][j] =
          a[0][j] * cinp[i][0] + a[1][j] * cinp[i][1] + a[2][j] * cinp[i][2];
    }
  }
}

void fphi_to_cphi(PMEUnit pme_u, const real (*fphi)[20], real (*cphi)[10]) {
  auto& st = pme_u.obj();
  int nfft1 = st.nfft1;
  int nfft2 = st.nfft2;
  int nfft3 = st.nfft3;

  // data qi1  / 1, 2, 3, 1, 1, 2 /
  // data qi2  / 1, 2, 3, 2, 3, 3 /
  constexpr int qi1[] = {0, 1, 2, 0, 0, 1};
  constexpr int qi2[] = {0, 1, 2, 1, 2, 2};

  real ftc3[3][3];
  real ftc6[6][6];
  real a[3][3];

  #pragma acc parallel loop deviceptr(box,fphi,cphi)\
              private(ftc3[0:3][0:3],ftc6[0:6][0:6],a[0:3][0:3])
  for (int iatom = 0; iatom < n; ++iatom) {

    // see also subroutine frac_to_cart in pmestuf.f
    // set the reciprocal vector transformation matrix

    a[0][0] = nfft1 * box->recip[0][0];
    a[1][0] = nfft2 * box->recip[1][0];
    a[2][0] = nfft3 * box->recip[2][0];
    a[0][1] = nfft1 * box->recip[0][1];
    a[1][1] = nfft2 * box->recip[1][1];
    a[2][1] = nfft3 * box->recip[2][1];
    a[0][2] = nfft1 * box->recip[0][2];
    a[1][2] = nfft2 * box->recip[1][2];
    a[2][2] = nfft3 * box->recip[2][2];

    // get the fractional to Cartesian conversion matrix

    #pragma acc loop independent collapse(2)
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        ftc3[j][i] = a[j][i];
      }
    }
    #pragma acc loop independent
    for (int i1 = 0; i1 < 3; ++i1) {
      int k = qi1[i1];
      #pragma acc loop independent
      for (int i2 = 0; i2 < 3; ++i2) {
        int i = qi1[i2];
        ftc6[i2][i1] = a[i][k] * a[i][k];
      }
      #pragma acc loop independent
      for (int i2 = 3; i2 < 6; ++i2) {
        int i = qi1[i2];
        int j = qi2[i2];
        ftc6[i2][i1] = 2 * a[i][k] * a[j][k];
      }
    }

    #pragma acc loop independent
    for (int i1 = 3; i1 < 6; ++i1) {
      int k = qi1[i1];
      int m = qi2[i1];
      #pragma acc loop independent
      for (int i2 = 0; i2 < 3; ++i2) {
        int i = qi1[i2];
        ftc6[i2][i1] = a[i][k] * a[i][m];
      }
      #pragma acc loop independent
      for (int i2 = 3; i2 < 6; ++i2) {
        int i = qi1[i2];
        int j = qi2[i2];
        ftc6[i2][i1] = a[i][k] * a[j][m] + a[i][m] * a[j][k];
      }
    }

    // apply the transformation to get the Cartesian potential

    cphi[iatom][0] = fphi[iatom][0];
    #pragma acc loop independent
    for (int j = 1; j < 10; ++j) {
      cphi[iatom][j] = 0;
    }
    #pragma acc loop seq collapse(2)
    for (int j = 1; j < 4; ++j) {
      for (int k = 1; k < 4; ++k) {
        cphi[iatom][j] += ftc3[k - 1][j - 1] * fphi[iatom][k];
      }
    }
    #pragma acc loop seq collapse(2)
    for (int j = 4; j < 10; ++j) {
      for (int k = 4; k < 10; ++k) {
        cphi[iatom][j] += ftc6[k - 4][j - 4] * fphi[iatom][k];
      }
    }
  }
}
TINKER_NAMESPACE_END
