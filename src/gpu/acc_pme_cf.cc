/**
 * @file
 * Conversion between Cartesian and fractional coordinates.
 */

#include "gpu/acc_fmat.h"
#include "gpu/decl_box.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_pme.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void cmp_to_fmp(const real (*_cmp)[10], real (*_fmp)[10], int pme_unit) {
  pme_st& st = pme_obj(pme_unit);
  int nfft1 = st.nfft1;
  int nfft2 = st.nfft2;
  int nfft3 = st.nfft3;

  constexpr int qi1[] = {1, 2, 3, 1, 1, 2};
  constexpr int qi2[] = {1, 2, 3, 2, 3, 3};

  #pragma acc data deviceptr(box,_cmp,_fmp)
  #pragma acc parallel loop
  for (int iatom = 0; iatom < n; ++iatom) {
    real ctf_cpp[10][10];
    freal_mat<10> ctf(ctf_cpp);

    // see also subroutine cart_to_frac in pmestuf.f
    {
      real a_cpp[3][3];
      freal_mat3 a(a_cpp);
      freal_mat3 recip(box->recip);

      // set the reciprocal vector transformation matrix

      #pragma acc loop independent
      for (int i = 1; i <= 3; ++i) {
        a(1, i) = nfft1 * recip(i, 1);
        a(2, i) = nfft2 * recip(i, 2);
        a(3, i) = nfft3 * recip(i, 3);
      }

      // get the Cartesian to fractional conversion matrix

      real* ctfptr = &ctf_cpp[0][0];
      #pragma acc loop independent
      for (int i = 0; i < 10 * 10; ++i) {
        ctfptr[i] = 0;
      }
      ctf(1, 1) = 1;
      #pragma acc loop independent collapse(2)
      for (int i = 2; i <= 4; ++i) {
        for (int j = 2; j <= 4; ++j) {
          ctf(i, j) = a(i - 1, j - 1);
        }
      }
      #pragma acc loop independent
      for (int i1 = 1; i1 <= 3; ++i1) {
        int k = qi1[i1 - 1];
        #pragma acc loop independent
        for (int i2 = 1; i2 <= 6; ++i2) {
          int i = qi1[i2 - 1];
          int j = qi2[i2 - 1];
          ctf(i1 + 4, i2 + 4) = a(k, i) * a(k, j);
        }
      }
      #pragma acc loop independent
      for (int i1 = 1; i1 <= 3; ++i1) {
        int k = qi1[i1 - 1];
        int m = qi2[i1 - 1];
        #pragma acc loop independent
        for (int i2 = 1; i2 <= 6; ++i2) {
          int i = qi1[i2 - 1];
          int j = qi2[i2 - 1];
          ctf(i1 + 4, i2 + 4) = a(k, i) * a(m, j) + a(k, j) * a(m, i);
        }
      }
    }

    // apply the transformation to get the fractional multipoles
    {
      freal_mat<10> cmp(_cmp);
      freal_mat<10> fmp(_fmp);
      const int i = iatom;
      fmp(1, i) = ctf(1, 1) * cmp(1, i);
      #pragma acc loop independent
      for (int j = 2; j <= 10; ++j) {
        fmp(j, i) = 0;
      }
      #pragma acc loop seq collapse(2)
      for (int j = 2; j <= 4; ++j) {
        for (int k = 2; k <= 4; ++k) {
          fmp(j, i) += ctf(j, k) * cmp(k, i);
        }
      }
      #pragma acc loop seq collapse(2)
      for (int j = 5; j <= 10; ++j) {
        for (int k = 5; k <= 10; ++k) {
          fmp(j, i) += ctf(j, k) * cmp(k, i);
        }
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
    freal_mat<10> ftc(ftc_cpp);

    // see also subroutine frac_to_cart in pmestuf.f
    {
      real a_cpp[3][3];
      freal_mat3 a(a_cpp);
      freal_mat3 recip(box->recip);

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
      freal_mat<20> fphi(_fphi);
      freal_mat<10> cphi(_cphi);
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
