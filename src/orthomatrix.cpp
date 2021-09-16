#include "tool/orthomatrix.h"
#include "mathfunc_sinhc.h"
#include <algorithm>

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


//====================================================================//


template <class T>
static T get_eps()
{
   return (T)1.0e-3;
}


/**
 * \f[
 * e_2(a,b)=\frac{\exp(a)-\exp(b)}{a-b}.
 * \f]
 */
template <class T>
static T e2(T a, T b)
{
   T eps = get_eps<T>();
   if (std::fabs(a - b) <= eps) {
      T x = (a + b) / 2, y = (a - b) / 2;
      return std::exp(x) * sinhc(y);
   } else {
      return (std::exp(a) - std::exp(b)) / (a - b);
   }
}


/**
 * \f[
 * f_2(a,b)=\frac{f_1(a)-f_1(b)}{a-b}.
 * \f]
 */
template <class T>
static T f2(T a, T b)
{
   T eps = get_eps<T>();
   T zero = 0;
   if (std::fabs(a - b) <= 2 * eps) {
      if (std::fabs(a) <= eps and std::fabs(b) <= eps) {
         T k0 = 1. / 2, k1 = 1. / 6, k2 = 1. / 24;
         T k3 = 1. / 120, k4 = 1. / 720, k5 = 1. / 5040;
         T l0, l1, l2, l3, l4, l5;
         T a2 = a * a, a3 = a2 * a;
         T b2 = b * b, b3 = b2 * b;
         T ab = a * b;
         l0 = 1;
         l1 = a + b;
         l2 = a2 + ab + b2;
         l3 = a3 + a2 * b + a * b2 + b3;
         l4 = a2 * a2 + a3 * b + a2 * b2 + a * b3 + b2 * b2;
         l5 = a3 * a2 + a2 * a2 * b + a3 * b2 + a2 * b3 + a * b2 * b2 + b2 * b3;
         return l0 * k0 + l1 * k1 + l2 * k2 + l3 * k3 + l4 * k4 + l5 * k5;
      } else {
         T x, y; // fabs(x) >= fabs(y)
         if (std::fabs(a) >= std::fabs(b))
            x = a, y = b;
         else
            x = b, y = a;
         return (e2(x, y) - e2(zero, y)) / x;
      }
   } else {
      return (e2(a, zero) - e2(b, zero)) / (a - b); // e2(a,0) = f1(a)
   }
}


/**
 * \f[
 * e_3(a,b,c)=\exp(a)f_2(b-a,c-a).
 * \f]
 */
template <class T>
static T e3(T a, T b, T c)
{
   T arr[] = {a, b, c};
   std::sort(arr, arr + 3);
   a = arr[2], b = arr[1], c = arr[0];
   return std::exp(a) * f2(b - a, c - a);
}


/**
 * \f[
 * f_3(a,b,c)=\frac{f_2(a,b)-f_2(b,c)}{a-c}.
 * \f]
 */
template <class T>
static T f3(T a, T b, T c)
{
   T arr[] = {a, b, c};
   std::sort(arr, arr + 3);
   a = arr[2], b = arr[1], c = arr[0]; // a >= b >= c
   T eps = get_eps<T>();
   if ((a - c) <= 2 * eps) {
      if (std::fabs(a) <= eps and std::fabs(b) <= eps and std::fabs(c) <= eps) {
         T a2 = a * a, a3 = a2 * a;
         T b2 = b * b, b3 = b2 * b;
         T c2 = c * c, c3 = c2 * c;
         T k0 = 1. / 6, k1 = 1. / 24;
         T k2 = 1. / 120, k3 = 1. / 720, k4 = 1. / 5040;
         T l0, l1, l2, l3, l4;
         l0 = 1;
         l1 = a + b + c;
         l2 = a2 + a * b + a * c + b2 + b * c + c2;
         l3 = a3 + a2 * b + a * b2 + a * b * c + a2 * c + a * c2 + b3 + b2 * c +
            b * c2 + c3;
         l4 = a2 * a2 + a3 * b + a2 * b2 + a * b3 + b2 * b2 + a3 * c +
            a2 * b * c + a * b2 * c + b3 * c + a2 * c2 + a * b * c2 + b2 * c2 +
            a * c3 + b * c3 + c2 * c2;
         return l0 * k0 + l1 * k1 + l2 * k2 + l3 * k3 + l4 * k4;
      } else {
         T x, y, z; // fabs(x) = max(fabs(a),fabs(b),fabs(c))
         T fa = std::fabs(a), fb = std::fabs(b), fc = std::fabs(c);
         if (fa >= fb and fa >= fc)
            x = a, y = b, z = c;
         else if (fb >= fc and fb >= fa)
            x = b, y = c, z = a;
         else
            x = c, y = a, z = b;
         return (e3(x, y, z) - f2(y, z)) / x; // f2(y,z) = e3(0,y,z)
      }
   } else {
      return (f2(a, b) - f2(b, c)) / (a - c);
   }
}


template <class T>
void trimat_exp(T ans[3][3], T m[3][3], T t)
{
   constexpr int unknown = 0x00;
   constexpr int ortho = 0x01;
   constexpr int mono = 0x02;
   constexpr int tri = 0x04;
   int matrixtype = unknown;
   if (m[1][0] == (T)0 and m[2][0] == (T)0 and m[2][1] == (T)0) {
      matrixtype += tri;
      if (m[0][1] == (T)0 and m[1][2] == (T)0) {
         matrixtype += mono;
         if (m[0][2] == (T)0) {
            matrixtype += ortho;
         }
      }
   }

   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         ans[i][j] = 0;
   if (matrixtype & (tri + mono + ortho)) {
      T x = m[0][0] * t, y = m[1][1] * t, z = m[2][2] * t;
      ans[0][0] += std::exp(x);
      ans[1][1] += std::exp(y);
      ans[2][2] += std::exp(z);

      if (matrixtype & (tri + mono)) {
         T w = m[0][2] * t;
         ans[0][2] += w * e2(z, x);

         if (matrixtype & tri) {
            T u = m[0][1] * t, v = m[1][2] * t;
            ans[0][1] += u * e2(x, y);
            ans[1][2] += v * e2(y, z);
            ans[0][2] += u * v * e3(x, y, z);
         }
      }
   }
}
template void trimat_exp(double ans[3][3], double m[3][3], double t);
template void trimat_exp(float ans[3][3], float m[3][3], float t);


template <class T>
void trimat_expm1c(T ans[3][3], T m[3][3], T t)
{
   constexpr int unknown = 0x00;
   constexpr int ortho = 0x01;
   constexpr int mono = 0x02;
   constexpr int tri = 0x04;
   int matrixtype = unknown;
   if (m[1][0] == (T)0 and m[2][0] == (T)0 and m[2][1] == (T)0) {
      matrixtype += tri;
      if (m[0][1] == (T)0 and m[1][2] == (T)0) {
         matrixtype += mono;
         if (m[0][2] == (T)0) {
            matrixtype += ortho;
         }
      }
   }

   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         ans[i][j] = 0;
   if (matrixtype & (tri + mono + ortho)) {
      T zero = 0;
      T x = m[0][0] * t, y = m[1][1] * t, z = m[2][2] * t;
      ans[0][0] += e2(x, zero); // f1(x)
      ans[1][1] += e2(y, zero);
      ans[2][2] += e2(z, zero);

      if (matrixtype & (tri + mono)) {
         T w = m[0][2] * t;
         ans[0][2] += w * f2(z, x);

         if (matrixtype & tri) {
            T u = m[0][1] * t, v = m[1][2] * t;
            ans[0][1] += u * f2(x, y);
            ans[1][2] += v * f2(y, z);
            ans[0][2] += u * v * f3(x, y, z);
         }
      }
   }
}
template void trimat_expm1c(double ans[3][3], double m[3][3], double t);
template void trimat_expm1c(float ans[3][3], float m[3][3], float t);
}
