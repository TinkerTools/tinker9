#pragma once
#include "macro.h"
#include "math/libfunc.h"

namespace tinker {
#pragma acc routine seq
template <int n, class T = double>
void ludcmp(T* lu, int* indx)
{
   const T TINY = 1.0e-16;
   T big, temp;
   T vv[n];

   for (int i = 0; i < n; i++) {
      big = 0;
      for (int j = 0; j < n; j++) {
         if ((temp = REAL_ABS(lu[i * n + j])) > big)
            big = temp;
      }
      // if (big == 0) throw("Singular matrix in LUdcmp");
      vv[i] = 1 / big;
   }

   for (int k = 0; k < n; k++) {
      big = 0;
      int imax = k;
      for (int i = k; i < n; i++) {
         temp = vv[i] * REAL_ABS(lu[i * n + k]);
         if (temp > big) {
            big = temp;
            imax = i;
         }
      }
      if (k != imax) {
         for (int j = 0; j < n; j++) {
            temp = lu[imax * n + j];
            lu[imax * n + j] = lu[k * n + j];
            lu[k * n + j] = temp;
         }
         vv[imax] = vv[k];
      }
      indx[k] = imax;
      if (lu[k * n + k] == 0)
         lu[k * n + k] = TINY;
      for (int i = k + 1; i < n; i++) {
         temp = lu[i * n + k] /= lu[k * n + k];
         for (int j = k + 1; j < n; j++)
            lu[i * n + j] -= temp * lu[k * n + j];
      }
   }
}

#pragma acc routine seq
template <int n, class T = double>
void lubksb(const T* lu, T* x, const int* indx)
{
   int ii = 0;
   T sum;
   for (int i = 0; i < n; i++) {
      int ip = indx[i];
      sum = x[ip];
      x[ip] = x[i];
      if (ii != 0) {
         for (int j = ii - 1; j < i; j++)
            sum -= lu[i * n + j] * x[j];
      } else if (sum != 0) {
         ii = i + 1;
      }
      x[i] = sum;
   }
   for (int i = n - 1; i >= 0; i--) {
      sum = x[i];
      for (int j = i + 1; j < n; j++)
         sum -= lu[i * n + j] * x[j];
      x[i] = sum / lu[i * n + i];
   }
}

#pragma acc routine seq
template <int n, class R, class T = double>
void lumprove(const T* lu, T* x, const int* indx, const T* aref, const R* b)
{
   T r[n];
   for (int i = 0; i < n; i++) {
      double sdp = -b[i];
      for (int j = 0; j < n; j++)
         sdp += (double)aref[i * n + j] * (double)x[j];
      r[i] = sdp;
   }
   lubksb<n>(lu, r, indx);
   for (int i = 0; i < n; i++)
      x[i] -= r[i];
}

/**
 * \ingroup math
 * \brief This subroutine uses the LU decomposition method to solve the linear
 * system Ax = b, returning x in b. A is an n by n real symmetric matrix with
 * its upper triangle (including the diagonal) stored by rows.
 *
 * Literature reference:
 *    - <a href="http://numerical.recipes">
 *    W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery,
 *    Numerical Recipes (C++), 3rd Ed., Section 2.3,
 *    Cambridge University Press (2007).
 *    </a>
 */
template <int n, class R>
void symlusolve(const R* aUpRowMajor, R* b)
{
   #pragma acc serial async deviceptr(aUpRowMajor,b)
   {
      // A x = b
      double a[n][n];  // full A matrix
      double lu[n][n]; // LU decompostion
      double x[n];     // solution
      int indx[n];     // permutation

      // Initialize the internal arrays.
      int c = 0;
      for (int k = 0; k < n; ++k) {
         x[k] = b[k];
         for (int m = 0; m < k; ++m) {
            a[k][m] = a[m][k];
            lu[k][m] = a[k][m];
         }
         for (int m = k; m < n; ++m) {
            a[k][m] = aUpRowMajor[c];
            lu[k][m] = a[k][m];
            ++c;
         }
      }

      double* plu = &lu[0][0];
      double* pa = &a[0][0];
      ludcmp<n>(plu, indx);
      lubksb<n>(plu, x, indx);
      lumprove<n>(plu, x, indx, pa, b);

      // Copy out the solution.
      for (int k = 0; k < n; ++k) {
         b[k] = x[k];
      }
   }
}
}
