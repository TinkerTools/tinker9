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
void matmul3(T S[3][3], T A[3][3], T B[3][3]);


/**
 * C-Style
 *          [[11, 12, 13],
 * matrix =  [21, 22, 23],
 *           [31, 32, 33]]
 *
 * answer[3][3] = A[3][3] answer[3][3]
 */
template <class T>
void matmul3(T answer[3][3], T A[3][3]);


struct SymmMatrix
{
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
   static int solve(T A0[3][3], T Q[3][3], T w[3]);


   // S = O DIAG O^T
   template <class T>
   static void ODOt(T s[3][3], T o[3][3], T d[3]);
};
}
