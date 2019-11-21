#include "box.h"
#include "e_mpole.h"
#include "md.h"
#include "pme.h"

/**
 * @file
 * Conversion between Cartesian and fractional coordinates.
 */

TINKER_NAMESPACE_BEGIN
void rpole_to_cmp()
{

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

void cmp_to_fmp(PMEUnit pme_u, const real (*cmp)[10], real (*fmp)[10])
{
   auto& st = *pme_u;
   int nfft1 = st.nfft1;
   int nfft2 = st.nfft2;
   int nfft3 = st.nfft3;

   #pragma acc parallel loop deviceptr(box,cmp,fmp)
   for (int iatom = 0; iatom < n; ++iatom) {

      real a[3][3];

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

      // data qi1  / 1, 2, 3, 1, 1, 2 /
      // data qi2  / 1, 2, 3, 2, 3, 3 /
      constexpr int qi1[] = {0, 1, 2, 0, 0, 1};
      constexpr int qi2[] = {0, 1, 2, 1, 2, 2};

      real ctf6[6][6];

      // get the Cartesian to fractional conversion matrix

      #pragma acc loop seq
      for (int i1 = 0; i1 < 3; ++i1) {
         int k = qi1[i1];
         #pragma acc loop seq
         for (int i2 = 0; i2 < 6; ++i2) {
            int i = qi1[i2];
            int j = qi2[i2];
            ctf6[i2][i1] = a[i][k] * a[j][k];
         }
      }
      #pragma acc loop seq
      for (int i1 = 3; i1 < 6; ++i1) {
         int k = qi1[i1];
         int m = qi2[i1];
         #pragma acc loop seq
         for (int i2 = 0; i2 < 6; ++i2) {
            int i = qi1[i2];
            int j = qi2[i2];
            ctf6[i2][i1] = a[i][k] * a[j][m] + a[j][k] * a[i][m];
         }
      }

      // apply the transformation to get the fractional multipoles

      real cmpi[10], fmpi[10];
      cmpi[0] = cmp[iatom][0];
      cmpi[1] = cmp[iatom][1];
      cmpi[2] = cmp[iatom][2];
      cmpi[3] = cmp[iatom][3];
      cmpi[4] = cmp[iatom][4];
      cmpi[5] = cmp[iatom][5];
      cmpi[6] = cmp[iatom][6];
      cmpi[7] = cmp[iatom][7];
      cmpi[8] = cmp[iatom][8];
      cmpi[9] = cmp[iatom][9];

      fmpi[0] = cmpi[0];
      #pragma acc loop seq
      for (int j = 1; j < 4; ++j) {
         fmpi[j] = 0;
         #pragma acc loop seq
         for (int k = 1; k < 4; ++k) {
            fmpi[j] += a[k - 1][j - 1] * cmpi[k];
         }
      }
      #pragma acc loop seq
      for (int j = 4; j < 10; ++j) {
         fmpi[j] = 0;
         #pragma acc loop seq
         for (int k = 4; k < 10; ++k) {
            fmpi[j] += ctf6[k - 4][j - 4] * cmpi[k];
         }
      }

      fmp[iatom][0] = fmpi[0];
      fmp[iatom][1] = fmpi[1];
      fmp[iatom][2] = fmpi[2];
      fmp[iatom][3] = fmpi[3];
      fmp[iatom][4] = fmpi[4];
      fmp[iatom][5] = fmpi[5];
      fmp[iatom][6] = fmpi[6];
      fmp[iatom][7] = fmpi[7];
      fmp[iatom][8] = fmpi[8];
      fmp[iatom][9] = fmpi[9];
   }
}

void cuind_to_fuind(PMEUnit pme_u, const real (*cind)[3], const real (*cinp)[3],
                    real (*fuind)[3], real (*fuinp)[3])
{
   auto& st = *pme_u;
   int nfft1 = st.nfft1;
   int nfft2 = st.nfft2;
   int nfft3 = st.nfft3;

   real a[3][3];

   #pragma acc parallel loop independent deviceptr(box,\
               cind,cinp,fuind,fuinp) private(a[0:3][0:3])
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

void fphi_to_cphi(PMEUnit pme_u, const real (*fphi)[20], real (*cphi)[10])
{
   auto& st = *pme_u;
   int nfft1 = st.nfft1;
   int nfft2 = st.nfft2;
   int nfft3 = st.nfft3;

   #pragma acc parallel loop deviceptr(box,fphi,cphi)
   for (int iatom = 0; iatom < n; ++iatom) {

      real a[3][3];

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

      // data qi1  / 1, 2, 3, 1, 1, 2 /
      // data qi2  / 1, 2, 3, 2, 3, 3 /
      constexpr int qi1[] = {0, 1, 2, 0, 0, 1};
      constexpr int qi2[] = {0, 1, 2, 1, 2, 2};

      real ftc6[6][6];

      // get the fractional to Cartesian conversion matrix

      #pragma acc loop seq
      for (int i1 = 0; i1 < 3; ++i1) {
         int k = qi1[i1];
         #pragma acc loop seq
         for (int i2 = 0; i2 < 3; ++i2) {
            int i = qi1[i2];
            ftc6[i2][i1] = a[i][k] * a[i][k];
         }
         #pragma acc loop seq
         for (int i2 = 3; i2 < 6; ++i2) {
            int i = qi1[i2];
            int j = qi2[i2];
            ftc6[i2][i1] = 2 * a[i][k] * a[j][k];
         }
      }

      #pragma acc loop seq
      for (int i1 = 3; i1 < 6; ++i1) {
         int k = qi1[i1];
         int m = qi2[i1];
         #pragma acc loop seq
         for (int i2 = 0; i2 < 3; ++i2) {
            int i = qi1[i2];
            ftc6[i2][i1] = a[i][k] * a[i][m];
         }
         #pragma acc loop seq
         for (int i2 = 3; i2 < 6; ++i2) {
            int i = qi1[i2];
            int j = qi2[i2];
            ftc6[i2][i1] = a[i][k] * a[j][m] + a[i][m] * a[j][k];
         }
      }

      // apply the transformation to get the Cartesian potential

      real fphii[10], cphii[10];
      fphii[0] = fphi[iatom][0];
      fphii[1] = fphi[iatom][1];
      fphii[2] = fphi[iatom][2];
      fphii[3] = fphi[iatom][3];
      fphii[4] = fphi[iatom][4];
      fphii[5] = fphi[iatom][5];
      fphii[6] = fphi[iatom][6];
      fphii[7] = fphi[iatom][7];
      fphii[8] = fphi[iatom][8];
      fphii[9] = fphi[iatom][9];

      cphii[0] = fphii[0];
      #pragma acc loop seq
      for (int j = 1; j < 4; ++j) {
         cphii[j] = 0;
         #pragma acc loop seq
         for (int k = 1; k < 4; ++k) {
            cphii[j] += a[k - 1][j - 1] * fphii[k];
         }
      }
      #pragma acc loop seq
      for (int j = 4; j < 10; ++j) {
         cphii[j] = 0;
         #pragma acc loop seq
         for (int k = 4; k < 10; ++k) {
            cphii[j] += ftc6[k - 4][j - 4] * fphii[k];
         }
      }

      cphi[iatom][0] = cphii[0];
      cphi[iatom][1] = cphii[1];
      cphi[iatom][2] = cphii[2];
      cphi[iatom][3] = cphii[3];
      cphi[iatom][4] = cphii[4];
      cphi[iatom][5] = cphii[5];
      cphi[iatom][6] = cphii[6];
      cphi[iatom][7] = cphii[7];
      cphi[iatom][8] = cphii[8];
      cphi[iatom][9] = cphii[9];
   }
}
TINKER_NAMESPACE_END
