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
void cmp_to_fmp(const real (*cmp)[10], real (*fmp)[10], int pme_unit) {
  pme_st& st = pme_obj(pme_unit);
  const pme_st* dptr = pme_deviceptr(pme_unit);

  int nfft1 = st.nfft1;
  int nfft2 = st.nfft2;
  int nfft3 = st.nfft3;

  constexpr int qi1[] = {1, 2, 3, 1, 1, 2};
  constexpr int qi2[] = {1, 2, 3, 2, 3, 3};

  #pragma acc data deviceptr(box,cmp,fmp,dptr)
  #pragma acc parallel loop
  for (int iatom = 0; iatom < n; ++iatom) {
    real ctf_cpp[10][10], a_cpp[3][3];
    freal<10, 10> ctf(ctf_cpp);
    freal33 a(a_cpp);
    freal33 recip(box->recip);

    // see also
    // subroutine cart_to_frac in pmestuf.f
    {
      // set the reciprocal vector transformation matrix

      #pragma acc loop independent
      for (int i = 1; i <= 3; ++i) {
        a(1, i) = nfft1 * recip(i, 1);
        a(2, i) = nfft2 * recip(i, 2);
        a(3, i) = nfft3 * recip(i, 3);
      }

      // get the Cartesian to fractional conversion matrix

      #pragma acc loop independent collapse(2)
      for (int i = 1; i <= 10; ++i) {
        for (int j = 1; j <= 10; ++j) {
          ctf(j, i) = 0;
        }
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
  }
}
}
TINKER_NAMESPACE_END
