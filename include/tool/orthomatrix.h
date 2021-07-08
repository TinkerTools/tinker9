#pragma once
#include <cmath>

namespace tinker {
/**
 * C-Style
 *          [[11, 12, 13],
 * matrix =  [21, 22, 23],
 *           [31, 32, 33]]
 *
 * S[3][3] = A[3][3] B[3][3]
 */
template <class T>
void matmul3(T S[3][3], T A[3][3], T B[3][3])
{
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         S[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
}


/**
 * C-Style
 *          [[11, 12, 13],
 * matrix =  [21, 22, 23],
 *           [31, 32, 33]]
 *
 * answer[3][3] = A[3][3] answer[3][3]
 */
template <class T>
void matmul3(T answer[3][3], T A[3][3])
{
   T B[3][3];
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         B[i][j] = answer[i][j];
   matmul3(answer, A, B);
}


struct SymmMatrix
{
private:
   template <class T>
   static T SQR(T x)
   {
      return x * x;
   }

public:
   /**
    * Q^T A Q = diag(w)
    * Q diag(w) Q^T = A
    *
    *      (0a)      (1a)      (2a)
    * v0 = (0b) v1 = (1b) v2 = (2b)
    *      (0c)      (1c)      (2c)
    *
    * Q v0 = w[0] v0
    * Q v1 = w[1] v1
    * Q v2 = w[2] v2
    *
    *           [[0a, 1a, 2a],
    * Q[3][3] =  [0b, 1b, 2b],
    *            [0c, 1c, 2c]]
    *
    * https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html
    *
    * Joachim Kopp
    * Efficient numerical diagonalization of hermitian 3x3 matrices
    * Int. J. Mod. Phys. C 19 (2008) 523-548
    * arXiv.org: physics/0610206
    */
   template <class T>
   static int solve(T A0[3][3], T Q[3][3], T w[3])
   {
      const int n = 3;
      T A[3][3];
      T sd, so;            // Sums of diagonal resp. off-diagonal elements
      T s, c;              // sin(phi), cos(phi), tan(phi)
      T t, g, h, z, theta; // Temporary storage
      T thresh;

      for (int i = 0; i < n; ++i)
         for (int j = 0; j < n; ++j)
            A[i][j] = A0[i][j];

      // Initialize Q to the identitity matrix
      for (int i = 0; i < n; i++) {
         Q[i][i] = 1.0;
         for (int j = 0; j < i; j++) {
            Q[i][j] = 0.0;
            Q[j][i] = 0.0;
         }
      }

      // Initialize w to diag(A)
      for (int i = 0; i < n; i++)
         w[i] = A[i][i];

      // Calculate SQR(tr(A))
      sd = 0.0;
      for (int i = 0; i < n; i++)
         sd += std::fabs(w[i]);
      sd = SQR(sd);

      // Main iteration loop
      for (int nIter = 0; nIter < 50; nIter++) {
         // Test for convergence
         so = 0.0;
         for (int p = 0; p < n; p++)
            for (int q = p + 1; q < n; q++)
               so += std::fabs(A[p][q]);
         if (so == 0.0)
            return 0;

         if (nIter < 4)
            thresh = 0.2 * so / SQR(n);
         else
            thresh = 0.0;

         // Do sweep
         for (int p = 0; p < n; p++)
            for (int q = p + 1; q < n; q++) {
               g = 100.0 * std::fabs(A[p][q]);
               if (nIter > 4 && std::fabs(w[p]) + g == std::fabs(w[p]) &&
                   std::fabs(w[q]) + g == std::fabs(w[q])) {
                  A[p][q] = 0.0;
               } else if (std::fabs(A[p][q]) > thresh) {
                  // Calculate Jacobi transformation
                  h = w[q] - w[p];
                  if (std::fabs(h) + g == std::fabs(h)) {
                     t = A[p][q] / h;
                  } else {
                     theta = 0.5 * h / A[p][q];
                     if (theta < 0.0)
                        t = -1.0 / (std::sqrt(1.0 + SQR(theta)) - theta);
                     else
                        t = 1.0 / (std::sqrt(1.0 + SQR(theta)) + theta);
                  }
                  c = 1.0 / std::sqrt(1.0 + SQR(t));
                  s = t * c;
                  z = t * A[p][q];

                  // Apply Jacobi transformation
                  A[p][q] = 0.0;
                  w[p] -= z;
                  w[q] += z;
                  for (int r = 0; r < p; r++) {
                     t = A[r][p];
                     A[r][p] = c * t - s * A[r][q];
                     A[r][q] = s * t + c * A[r][q];
                  }
                  for (int r = p + 1; r < q; r++) {
                     t = A[p][r];
                     A[p][r] = c * t - s * A[r][q];
                     A[r][q] = s * t + c * A[r][q];
                  }
                  for (int r = q + 1; r < n; r++) {
                     t = A[p][r];
                     A[p][r] = c * t - s * A[q][r];
                     A[q][r] = s * t + c * A[q][r];
                  }

                  // Update eigenvectors
                  for (int r = 0; r < n; r++) {
                     t = Q[r][p];
                     Q[r][p] = c * t - s * Q[r][q];
                     Q[r][q] = s * t + c * Q[r][q];
                  }
               }
            }
      }

      return -1;
   }
};
}
