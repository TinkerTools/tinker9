/**
 * @file
 * Conversion between Cartesian and fractional coordinates.
 */

#include "gpu/acc_fmat.h"
#include "gpu/decl_box.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_pme.h"
#include "gpu/e_mpole.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void cmp_to_fmp(real (*_fmp)[10], int pme_unit) {
  pme_st& st = pme_obj(pme_unit);
  int nfft1 = st.nfft1;
  int nfft2 = st.nfft2;
  int nfft3 = st.nfft3;

  constexpr int qi1[] = {1, 2, 3, 1, 1, 2};
  constexpr int qi2[] = {1, 2, 3, 2, 3, 3};

  #pragma acc parallel loop deviceptr(box,rpole,_fmp)
  for (int iatom = 0; iatom < n; ++iatom) {
    real ctf3_cpp[3][3];
    real ctf6_cpp[6][6];
    real a_cpp[3][3];

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

    real cmp_cpp[10] = {
        rpole[iatom][mpl_pme_0],      rpole[iatom][mpl_pme_x],
        rpole[iatom][mpl_pme_y],      rpole[iatom][mpl_pme_z],
        rpole[iatom][mpl_pme_xx],     rpole[iatom][mpl_pme_yy],
        rpole[iatom][mpl_pme_zz],     2 * rpole[iatom][mpl_pme_xy],
        2 * rpole[iatom][mpl_pme_xz], 2 * rpole[iatom][mpl_pme_yz]};
    farray_real cmp(cmp_cpp);
    fmat_real<10> fmp(_fmp);
    // to use the fortran matrix wrapper
    // change the 0-based atom number to 1-based
    const int i = iatom + 1;

    fmp(1, i) = cmp(1);
    #pragma acc loop independent
    for (int j = 2; j <= 10; ++j) {
      fmp(j, i) = 0;
    }
    #pragma acc loop seq collapse(2)
    for (int j = 2; j <= 4; ++j) {
      for (int k = 2; k <= 4; ++k) {
        fmp(j, i) += ctf3(j - 1, k - 1) * cmp(k);
      }
    }
    #pragma acc loop seq collapse(2)
    for (int j = 5; j <= 10; ++j) {
      for (int k = 5; k <= 10; ++k) {
        fmp(j, i) += ctf6(j - 4, k - 4) * cmp(k);
      }
    }
  }
}

void fphi_to_cphi(const real (*_fphi)[20], real (*_cphi)[10], int pme_unit) {
  pme_st& st = pme_obj(pme_unit);
  int nfft1 = st.nfft1;
  int nfft2 = st.nfft2;
  int nfft3 = st.nfft3;

  constexpr int qi1[] = {1, 2, 3, 1, 1, 2};
  constexpr int qi2[] = {1, 2, 3, 2, 3, 3};

  #pragma acc data deviceptr(box,_fphi,_cphi)
  #pragma acc parallel loop
  for (int iatom = 0; iatom < n; ++iatom) {
    real ftc_cpp[10][10];
    fmat_real<10> ftc(ftc_cpp);

    // see also subroutine frac_to_cart in pmestuf.f
    {
      real a_cpp[3][3];
      fmat_real3 a(a_cpp);
      fmat_real3 recip(box->recip);

      // set the reciprocal vector transformation matrix

      #pragma acc loop independent
      for (int i = 1; i <= 3; ++i) {
        a(i, 1) = nfft1 * recip(i, 1);
        a(i, 2) = nfft2 * recip(i, 2);
        a(i, 3) = nfft3 * recip(i, 3);
      }

      // get the fractional to Cartesian conversion matrix

      real* ftcptr = &ftc_cpp[0][0];
      #pragma acc loop independent
      for (int i = 0; i < 10 * 10; ++i) {
        ftcptr[i] = 0;
      }
      ftc(1, 1) = 1;
      #pragma acc loop independent collapse(2)
      for (int i = 2; i <= 4; ++i) {
        for (int j = 2; j <= 4; ++j) {
          ftc(i, j) = a(i - 1, j - 1);
        }
      }
      for (int i1 = 1; i1 <= 3; ++i1) {
        int k = qi1[i1 - 1];
        for (int i2 = 1; i2 <= 3; ++i2) {
          int i = qi1[i2 - 1];
          ftc(i1 + 4, i2 + 4) = a(k, i) * a(k, i);
        }
        for (int i2 = 4; i2 <= 6; ++i2) {
          int i = qi1[i2 - 1];
          int j = qi2[i2 - 1];
          ftc(i1 + 4, i2 + 4) = 2 * a(k, i) * a(k, j);
        }
      }
      for (int i1 = 4; i1 <= 6; ++i1) {
        int k = qi1[i1 - 1];
        int m = qi2[i1 - 1];
        for (int i2 = 1; i2 <= 3; ++i2) {
          int i = qi1[i2 - 1];
          ftc(i1 + 4, i2 + 4) = a(k, i) * a(m, i);
        }
        for (int i2 = 4; i2 <= 6; ++i2) {
          int i = qi1[i2 - 1];
          int j = qi2[i2 - 1];
          ftc(i1 + 4, i2 + 4) = a(k, i) * a(m, j) + a(m, i) * a(k, j);
        }
      }
    }

    // apply the transformation to get the Cartesian potential
    {
      fmat_real<20> fphi(_fphi);
      fmat_real<10> cphi(_cphi);
      const int i = iatom;
      cphi(1, i) = ftc(1, 1) * fphi(1, i);
      #pragma acc loop independent
      for (int j = 2; j <= 10; ++j) {
        cphi(j, i) = 0;
      }
      #pragma acc loop seq collapse(2)
      for (int j = 2; j <= 4; ++j) {
        for (int k = 2; k <= 4; ++k) {
          cphi(j, i) += ftc(j, k) * fphi(k, i);
        }
      }
      #pragma acc loop seq collapse(2)
      for (int j = 5; j <= 10; ++j) {
        for (int k = 5; k <= 10; ++k) {
          cphi(j, i) += ftc(j, k) * fphi(k, i);
        }
      }
    }
  }
}
}
TINKER_NAMESPACE_END
