#include "tool/orthomatrix.h"


namespace tinker {
template <class T>
void matmul3(T S[3][3], T A[3][3], T B[3][3])
{
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         S[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
}
template void matmul3(double s[3][3], double A[3][3], double B[3][3]);
template void matmul3(float s[3][3], float A[3][3], float B[3][3]);


template <class T>
void matmul3(T answer[3][3], T A[3][3])
{
   T B[3][3];
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         B[i][j] = answer[i][j];
   matmul3(answer, A, B);
}
template void matmul3(double ans[3][3], double A[3][3]);
template void matmul3(float ans[3][3], float A[3][3]);


template <class T>
static T SQR(T x)
{
   return x * x;
}
template <class T>
int SymmMatrix::solve(T A0[3][3], T Q[3][3], T w[3])
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
template int SymmMatrix::solve(double A0[3][3], double Q[3][3], double w[3]);
template int SymmMatrix::solve(float A0[3][3], float Q[3][3], float w[3]);


template <class T>
void SymmMatrix::ODOt(T s[3][3], T o[3][3], T d[3])
{
   // clang-format off
   s[0][0] = d[0]*o[0][0]*o[0][0] + d[1]*o[0][1]*o[0][1] + d[2]*o[0][2]*o[0][2];
   s[0][1] = d[0]*o[0][0]*o[1][0] + d[1]*o[0][1]*o[1][1] + d[2]*o[0][2]*o[1][2];
   s[0][2] = d[0]*o[0][0]*o[2][0] + d[1]*o[0][1]*o[2][1] + d[2]*o[0][2]*o[2][2];
   // s[1][0] = d[0]*o[0][0]*o[1][0] + d[1]*o[0][1]*o[1][1] + d[2]*o[0][2]*o[1][2];
   s[1][0] = s[0][1];
   s[1][1] = d[0]*o[1][0]*o[1][0] + d[1]*o[1][1]*o[1][1] + d[2]*o[1][2]*o[1][2];
   s[1][2] = d[0]*o[1][0]*o[2][0] + d[1]*o[1][1]*o[2][1] + d[2]*o[1][2]*o[2][2];
   // s[2][0] = d[0]*o[0][0]*o[2][0] + d[1]*o[0][1]*o[2][1] + d[2]*o[0][2]*o[2][2];
   // s[2][1] = d[0]*o[1][0]*o[2][0] + d[1]*o[1][1]*o[2][1] + d[2]*o[1][2]*o[2][2];
   s[2][0] = s[0][2];
   s[2][1] = s[1][2];
   s[2][2] = d[0]*o[2][0]*o[2][0] + d[1]*o[2][1]*o[2][1] + d[2]*o[2][2]*o[2][2];
   // clang-format on
}
template void SymmMatrix::ODOt(double s[3][3], double o[3][3], double d[3]);
template void SymmMatrix::ODOt(float s[3][3], float o[3][3], float d[3]);
}
