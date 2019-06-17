#include "gpu/acc_fmat.h"
#include "gpu/decl_box.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_pme.h"
#include "gpu/e_mpole.h"

/**
 * @file
 * Conversion between Cartesian and fractional coordinates.
 */

TINKER_NAMESPACE_BEGIN
namespace gpu {
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

void cmp_to_fmp(int pme_unit, const real (*cmp)[10], real (*_fmp)[10]) {
  pme_st& st = pme_obj(pme_unit);
  int nfft1 = st.nfft1;
  int nfft2 = st.nfft2;
  int nfft3 = st.nfft3;

  constexpr int qi1[] = {1, 2, 3, 1, 1, 2};
  constexpr int qi2[] = {1, 2, 3, 2, 3, 3};

  real ctf3_cpp[3][3];
  real ctf6_cpp[6][6];
  real a_cpp[3][3];

  #pragma acc parallel loop deviceptr(box,cmp,_fmp)\
              private(ctf3_cpp[0:3][0:3],ctf6_cpp[0:6][0:6],a_cpp[0:3][0:3])
  for (int iatom = 0; iatom < n; ++iatom) {

    // see also subroutine cart_to_frac in pmestuf.f

    fmat_real<3> ctf3(ctf3_cpp);
    fmat_real<6> ctf6(ctf6_cpp);
    fmat_real3 a(a_cpp);
    fmat_real3 recip(box->recip);

    // set the reciprocal vector transformation matrix

    a(1, 1) = nfft1 * recip(1, 1);
    a(2, 1) = nfft2 * recip(1, 2);
    a(3, 1) = nfft3 * recip(1, 3);
    a(1, 2) = nfft1 * recip(2, 1);
    a(2, 2) = nfft2 * recip(2, 2);
    a(3, 2) = nfft3 * recip(2, 3);
    a(1, 3) = nfft1 * recip(3, 1);
    a(2, 3) = nfft2 * recip(3, 2);
    a(3, 3) = nfft3 * recip(3, 3);

    // get the Cartesian to fractional conversion matrix

    #pragma acc loop independent collapse(2)
    for (int i = 1; i <= 3; ++i) {
      for (int j = 1; j <= 3; ++j) {
        ctf3(i, j) = a(i, j);
      }
    }
    #pragma acc loop independent
    for (int i1 = 1; i1 <= 3; ++i1) {
      int k = qi1[i1 - 1];
      #pragma acc loop independent
      for (int i2 = 1; i2 <= 6; ++i2) {
        int i = qi1[i2 - 1];
        int j = qi2[i2 - 1];
        ctf6(i1, i2) = a(k, i) * a(k, j);
      }
    }
    #pragma acc loop independent
    for (int i1 = 4; i1 <= 6; ++i1) {
      int k = qi1[i1 - 1];
      int m = qi2[i1 - 1];
      #pragma acc loop independent
      for (int i2 = 1; i2 <= 6; ++i2) {
        int i = qi1[i2 - 1];
        int j = qi2[i2 - 1];
        ctf6(i1, i2) = a(k, i) * a(m, j) + a(k, j) * a(m, i);
      }
    }

    // apply the transformation to get the fractional multipoles

    fmat_real<10> fmp(_fmp);

    // to use the fortran matrix wrapper
    // change the 0-based atom number to 1-based
    const int i = iatom + 1;

    fmp(1, i) = cmp[iatom][0];
    #pragma acc loop independent
    for (int j = 2; j <= 10; ++j) {
      fmp(j, i) = 0;
    }
    #pragma acc loop seq collapse(2)
    for (int j = 2; j <= 4; ++j) {
      for (int k = 2; k <= 4; ++k) {
        fmp(j, i) += ctf3(j - 1, k - 1) * cmp[iatom][k - 1];
      }
    }
    #pragma acc loop seq collapse(2)
    for (int j = 5; j <= 10; ++j) {
      for (int k = 5; k <= 10; ++k) {
        fmp(j, i) += ctf6(j - 4, k - 4) * cmp[iatom][k - 1];
      }
    }
  }
}

void cuind_to_fuind(int pme_unit, const real (*cind)[3], const real (*cinp)[3],
                    real (*fuind)[3], real (*fuinp)[3]) {
  pme_st& st = pme_obj(pme_unit);
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

void fphi_to_cphi(int pme_unit, const real (*_fphi)[20], real (*_cphi)[10]) {
  pme_st& st = pme_obj(pme_unit);
  int nfft1 = st.nfft1;
  int nfft2 = st.nfft2;
  int nfft3 = st.nfft3;

  constexpr int qi1[] = {1, 2, 3, 1, 1, 2};
  constexpr int qi2[] = {1, 2, 3, 2, 3, 3};

  real ftc3_cpp[3][3];
  real ftc6_cpp[6][6];
  real a_cpp[3][3];

  #pragma acc parallel loop deviceptr(box,_fphi,_cphi)\
              private(ftc3_cpp[0:3][0:3],ftc6_cpp[0:6][0:6],a_cpp[0:3][0:3])
  for (int iatom = 0; iatom < n; ++iatom) {

    // see also subroutine frac_to_cart in pmestuf.f

    fmat_real<3> ftc3(ftc3_cpp);
    fmat_real<6> ftc6(ftc6_cpp);
    fmat_real3 a(a_cpp);
    fmat_real3 recip(box->recip);

    // set the reciprocal vector transformation matrix

    a(1, 1) = nfft1 * recip(1, 1);
    a(1, 2) = nfft2 * recip(1, 2);
    a(1, 3) = nfft3 * recip(1, 3);
    a(2, 1) = nfft1 * recip(2, 1);
    a(2, 2) = nfft2 * recip(2, 2);
    a(2, 3) = nfft3 * recip(2, 3);
    a(3, 1) = nfft1 * recip(3, 1);
    a(3, 2) = nfft2 * recip(3, 2);
    a(3, 3) = nfft3 * recip(3, 3);

    // get the fractional to Cartesian conversion matrix

    #pragma acc loop independent collapse(2)
    for (int i = 1; i <= 3; ++i) {
      for (int j = 1; j <= 3; ++j) {
        ftc3(i, j) = a(i, j);
      }
    }
    #pragma acc loop independent
    for (int i1 = 1; i1 <= 3; ++i1) {
      int k = qi1[i1 - 1];
      #pragma acc loop independent
      for (int i2 = 1; i2 <= 3; ++i2) {
        int i = qi1[i2 - 1];
        ftc6(i1, i2) = a(k, i) * a(k, i);
      }
      #pragma acc loop independent
      for (int i2 = 4; i2 <= 6; ++i2) {
        int i = qi1[i2 - 1];
        int j = qi2[i2 - 1];
        ftc6(i1, i2) = 2 * a(k, i) * a(k, j);
      }
    }
    #pragma acc loop independent
    for (int i1 = 4; i1 <= 6; ++i1) {
      int k = qi1[i1 - 1];
      int m = qi2[i1 - 1];
      #pragma acc loop independent
      for (int i2 = 1; i2 <= 3; ++i2) {
        int i = qi1[i2 - 1];
        ftc6(i1, i2) = a(k, i) * a(m, i);
      }
      #pragma acc loop independent
      for (int i2 = 4; i2 <= 6; ++i2) {
        int i = qi1[i2 - 1];
        int j = qi2[i2 - 1];
        ftc6(i1, i2) = a(k, i) * a(m, j) + a(m, i) * a(k, j);
      }
    }

    // apply the transformation to get the Cartesian potential

    fmat_real<20> fphi(_fphi);
    fmat_real<10> cphi(_cphi);
    // to use the fortran matrix wrapper
    // change the 0-based atom number to 1-based
    const int i = iatom + 1;

    cphi(1, i) = fphi(1, i);
    #pragma acc loop independent
    for (int j = 2; j <= 10; ++j) {
      cphi(j, i) = 0;
    }
    #pragma acc loop seq collapse(2)
    for (int j = 2; j <= 4; ++j) {
      for (int k = 2; k <= 4; ++k) {
        cphi(j, i) += ftc3(j - 1, k - 1) * fphi(k, i);
      }
    }
    #pragma acc loop seq collapse(2)
    for (int j = 5; j <= 10; ++j) {
      for (int k = 5; k <= 10; ++k) {
        cphi(j, i) += ftc6(j - 4, k - 4) * fphi(k, i);
      }
    }
  }
}
}
TINKER_NAMESPACE_END
