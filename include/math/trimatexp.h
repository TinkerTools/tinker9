#pragma once

namespace tinker {
/// \ingroup math
/// \brief Matrix multiplication of two 3 by 3 matrices. \f$ R = A B \f$.
///
/// Matrices are stored in the row-major order (C-style).
///
/// \param[out] R  Matrix for the result.
/// \param[in]  A  Matrix on the left.
/// \param[in]  B  Matrix on the right.
template <class T>
void matmul3(T R[3][3], T A[3][3], T B[3][3]);

/// \ingroup math
/// \brief Matrix multiplication of two 3 by 3 matrices. \f$ R = A R \f$.
///
/// Matrices are stored in the row-major order (C-style).
///
/// \param[in,out] R  Matrix to be multiplied on the left.
/// \param[in]     A  Matrix on the left.
template <class T>
void matmul3(T R[3][3], T A[3][3]);

/// \ingroup math
/// \brief \f$ \exp(mt) \f$. Matrix m is 3 by 3 upper triangular.
///
/// Matrices are stored in the row-major order (C-style).
///
/// \param[out] ans  3 by 3 matrix for the result.
/// \param[in]  m    3 by 3 upper triangular matrix.
/// \param[in]  t    Scalar to scale the matrix m.
///
/// \note If the input matrix is not an upper triangular matrix, the result is undefined.
template <class T>
void trimatExp(T ans[3][3], T m[3][3], T t);

/// \ingroup math
/// \brief \f$ (\exp(mt)-I)/(mt) \f$. Matrix m is 3 by 3 upper triangular.
///
/// Matrices are stored in the row-major order (C-style).
///
/// \param[out] ans  3 by 3 matrix for the result.
/// \param[in]  m    3 by 3 upper triangular matrix.
/// \param[in]  t    Scalar to scale the matrix m.
///
/// \note If the input matrix is not an upper triangular matrix, the result is undefined.
template <class T>
void trimatExpm1c(T ans[3][3], T m[3][3], T t);

/// \ingroup math
/// \brief \f$ t (\exp(mt)-I)/(mt) \f$. Matrix m is 3 by 3 upper triangular.
///
/// Matrices are stored in the row-major order (C-style).
///
/// \param[out] ans  3 by 3 matrix for the result.
/// \param[in]  m    3 by 3 upper triangular matrix.
/// \param[in]  t    Scalar to scale the matrix m.
///
/// \note If the input matrix is not an upper triangular matrix, the result is undefined.
template <class T>
void trimatTExpm1c(T ans[3][3], T m[3][3], T t)
{
   trimatExpm1c(ans, m, t);
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         ans[i][j] *= t;
}
}
