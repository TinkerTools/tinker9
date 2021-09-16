#include "test.h"
#include "tool/orthomatrix.h"
using namespace tinker;

TEST_CASE("OrthoMatrix", "[util][math][orthomatrix]")
{
   SECTION("Example 1")
   {
      double A[3][3] = {{2.0, 2.0 / 3.0, -2.0 / 3.0},
                        {2.0 / 3.0, 5.0 / 3.0, 0.0},
                        {-2.0 / 3.0, 0.0, 7.0 / 3.0}},
             o[3][3] = {0.0}, w[3] = {0.0};
      SymmMatrix::solve(A, o, w);

      double diag[3][3] = {0}, ot[3][3];
      diag[0][0] = w[0];
      diag[1][1] = w[1];
      diag[2][2] = w[2];
      for (int i = 0; i < 3; ++i)
         for (int j = 0; j < 3; ++j)
            ot[i][j] = o[j][i];

      double tmp[3][3], a1[3][3], a2[3][3];
      double eps = 1.0e-8;
      // A = O diag O^T
      matmul3(tmp, diag, ot);
      matmul3(a1, o, tmp);
      REQUIRE(A[0][0] == Approx(a1[0][0]).margin(eps));
      REQUIRE(A[0][1] == Approx(a1[0][1]).margin(eps));
      REQUIRE(A[0][2] == Approx(a1[0][2]).margin(eps));
      REQUIRE(A[1][0] == Approx(a1[1][0]).margin(eps));
      REQUIRE(A[1][1] == Approx(a1[1][1]).margin(eps));
      REQUIRE(A[1][2] == Approx(a1[1][2]).margin(eps));
      REQUIRE(A[2][0] == Approx(a1[2][0]).margin(eps));
      REQUIRE(A[2][1] == Approx(a1[2][1]).margin(eps));
      REQUIRE(A[2][2] == Approx(a1[2][2]).margin(eps));

      // diag = O^T A O
      matmul3(tmp, ot, A);
      matmul3(a2, tmp, o);
      REQUIRE(diag[0][0] == Approx(a2[0][0]).margin(eps));
      REQUIRE(diag[0][1] == Approx(a2[0][1]).margin(eps));
      REQUIRE(diag[0][2] == Approx(a2[0][2]).margin(eps));
      REQUIRE(diag[1][0] == Approx(a2[1][0]).margin(eps));
      REQUIRE(diag[1][1] == Approx(a2[1][1]).margin(eps));
      REQUIRE(diag[1][2] == Approx(a2[1][2]).margin(eps));
      REQUIRE(diag[2][0] == Approx(a2[2][0]).margin(eps));
      REQUIRE(diag[2][1] == Approx(a2[2][1]).margin(eps));
      REQUIRE(diag[2][2] == Approx(a2[2][2]).margin(eps));
   }

   SECTION("Example 2")
   {
      double A[3][3] = {{19., 20., -16.}, {20., 13., 4.}, {-16., 4., 31.}},
             o[3][3] = {0.0}, w[3] = {0.0};
      SymmMatrix::solve(A, o, w);

      double diag[3][3] = {0}, ot[3][3];
      diag[0][0] = w[0];
      diag[1][1] = w[1];
      diag[2][2] = w[2];
      for (int i = 0; i < 3; ++i)
         for (int j = 0; j < 3; ++j)
            ot[i][j] = o[j][i];

      double tmp[3][3], a1[3][3], a2[3][3];
      double eps = 1.0e-8;
      // A = O diag O^T
      matmul3(tmp, diag, ot);
      matmul3(a1, o, tmp);
      REQUIRE(A[0][0] == Approx(a1[0][0]).margin(eps));
      REQUIRE(A[0][1] == Approx(a1[0][1]).margin(eps));
      REQUIRE(A[0][2] == Approx(a1[0][2]).margin(eps));
      REQUIRE(A[1][0] == Approx(a1[1][0]).margin(eps));
      REQUIRE(A[1][1] == Approx(a1[1][1]).margin(eps));
      REQUIRE(A[1][2] == Approx(a1[1][2]).margin(eps));
      REQUIRE(A[2][0] == Approx(a1[2][0]).margin(eps));
      REQUIRE(A[2][1] == Approx(a1[2][1]).margin(eps));
      REQUIRE(A[2][2] == Approx(a1[2][2]).margin(eps));

      // diag = O^T A O
      matmul3(tmp, ot, A);
      matmul3(a2, tmp, o);
      REQUIRE(diag[0][0] == Approx(a2[0][0]).margin(eps));
      REQUIRE(diag[0][1] == Approx(a2[0][1]).margin(eps));
      REQUIRE(diag[0][2] == Approx(a2[0][2]).margin(eps));
      REQUIRE(diag[1][0] == Approx(a2[1][0]).margin(eps));
      REQUIRE(diag[1][1] == Approx(a2[1][1]).margin(eps));
      REQUIRE(diag[1][2] == Approx(a2[1][2]).margin(eps));
      REQUIRE(diag[2][0] == Approx(a2[2][0]).margin(eps));
      REQUIRE(diag[2][1] == Approx(a2[2][1]).margin(eps));
      REQUIRE(diag[2][2] == Approx(a2[2][2]).margin(eps));
   }

   SECTION("ODOt")
   {
      double O[3][3] = {{2. / 3, -2. / 3, 1. / 3},
                        {1. / 3, 2. / 3, 2. / 3},
                        {2. / 3, 1. / 3, -2. / 3}};
      double Ot[3][3];
      for (int i = 0; i < 3; ++i)
         for (int j = 0; j < 3; ++j)
            Ot[i][j] = O[j][i];
      double diag[3] = {9., 18., 27.};
      double s[3][3], eps = 1.0e-8;

      double ref1[3][3] = {{15., 0., -6.}, {0., 21., -6.}, {-6., -6., 18.}};
      SymmMatrix::ODOt(s, O, diag);
      REQUIRE(s[0][0] == Approx(ref1[0][0]).margin(eps));
      REQUIRE(s[0][1] == Approx(ref1[0][1]).margin(eps));
      REQUIRE(s[0][2] == Approx(ref1[0][2]).margin(eps));
      REQUIRE(s[1][0] == Approx(ref1[1][0]).margin(eps));
      REQUIRE(s[1][1] == Approx(ref1[1][1]).margin(eps));
      REQUIRE(s[1][2] == Approx(ref1[1][2]).margin(eps));
      REQUIRE(s[2][0] == Approx(ref1[2][0]).margin(eps));
      REQUIRE(s[2][1] == Approx(ref1[2][1]).margin(eps));
      REQUIRE(s[2][2] == Approx(ref1[2][2]).margin(eps));

      double ref2[3][3] = {{18., 6., -6.}, {6., 15., 0.}, {-6., 0., 21.}};
      SymmMatrix::ODOt(s, Ot, diag);
      REQUIRE(s[0][0] == Approx(ref2[0][0]).margin(eps));
      REQUIRE(s[0][1] == Approx(ref2[0][1]).margin(eps));
      REQUIRE(s[0][2] == Approx(ref2[0][2]).margin(eps));
      REQUIRE(s[1][0] == Approx(ref2[1][0]).margin(eps));
      REQUIRE(s[1][1] == Approx(ref2[1][1]).margin(eps));
      REQUIRE(s[1][2] == Approx(ref2[1][2]).margin(eps));
      REQUIRE(s[2][0] == Approx(ref2[2][0]).margin(eps));
      REQUIRE(s[2][1] == Approx(ref2[2][1]).margin(eps));
      REQUIRE(s[2][2] == Approx(ref2[2][2]).margin(eps));
   }

   SECTION("a")
   {
      double A[3][3] = {{3, 0, 8}, {0, 4, 0}, {8, 0, 5}}, o[3][3], w[3];
      SymmMatrix::solve(A, o, w);
      double s[3][3];
      SymmMatrix::ODOt(s, o, w);

      double B[3][3] = {{3, 0, 8}, {0, 4, 0}, {0, 0, 5}};
      matmul3(B, s);
   }
}


TEST_CASE("MatrixExp", "[util][math][matrixexp]")
{
   unsigned long long c = 0x080000000ull;

   SECTION("orthogonal")
   {
      double m[3][3] = {{1, 0, 0}, //
                        {0, 2, 0}, //
                        {0, 0, 1}};
      double r[3][3] = {{1.105170918, 0, 0}, //
                        {0, 1.221402758, 0}, //
                        {0, 0, 1.105170918}};

      double a[3][3] = {0};
      trimat_exp(a, m, 0.1);

      REQUIRE((long long)(c * r[0][0]) == (long long)(c * a[0][0]));
      REQUIRE((long long)(c * r[0][1]) == (long long)(c * a[0][1]));
      REQUIRE((long long)(c * r[0][2]) == (long long)(c * a[0][2]));
      REQUIRE((long long)(c * r[1][0]) == (long long)(c * a[1][0]));
      REQUIRE((long long)(c * r[1][1]) == (long long)(c * a[1][1]));
      REQUIRE((long long)(c * r[1][2]) == (long long)(c * a[1][2]));
      REQUIRE((long long)(c * r[2][0]) == (long long)(c * a[2][0]));
      REQUIRE((long long)(c * r[2][1]) == (long long)(c * a[2][1]));
      REQUIRE((long long)(c * r[2][2]) == (long long)(c * a[2][2]));
   }

   SECTION("monoclinic")
   {
      double m[3][3] = {{1, 0, 0.5}, //
                        {0, 2, 0},   //
                        {0, 0, 1}};
      double r[3][3] = {{1.105170918, 0, 0.05525854590}, //
                        {0, 1.221402758, 0},             //
                        {0, 0, 1.105170918}};

      double a[3][3] = {0};
      trimat_exp(a, m, 0.1);

      REQUIRE((long long)(c * r[0][0]) == (long long)(c * a[0][0]));
      REQUIRE((long long)(c * r[0][1]) == (long long)(c * a[0][1]));
      REQUIRE((long long)(c * r[0][2]) == (long long)(c * a[0][2]));
      REQUIRE((long long)(c * r[1][0]) == (long long)(c * a[1][0]));
      REQUIRE((long long)(c * r[1][1]) == (long long)(c * a[1][1]));
      REQUIRE((long long)(c * r[1][2]) == (long long)(c * a[1][2]));
      REQUIRE((long long)(c * r[2][0]) == (long long)(c * a[2][0]));
      REQUIRE((long long)(c * r[2][1]) == (long long)(c * a[2][1]));
      REQUIRE((long long)(c * r[2][2]) == (long long)(c * a[2][2]));
   }

   SECTION("triclinic1")
   {
      double m[3][3] = {{1, 1. / 3., 0.5}, //
                        {0, 1, 0.2},       //
                        {0, 0, 1}};
      double r[3][3] = {{1.105170918, 0.03683903060, 0.05562693621}, //
                        {0, 1.105170918, 0.02210341836},             //
                        {0, 0, 1.105170918}};

      double a[3][3] = {0};
      trimat_exp(a, m, 0.1);

      REQUIRE((long long)(c * r[0][0]) == (long long)(c * a[0][0]));
      REQUIRE((long long)(c * r[0][1]) == (long long)(c * a[0][1]));
      REQUIRE((long long)(c * r[0][2]) == (long long)(c * a[0][2]));
      REQUIRE((long long)(c * r[1][0]) == (long long)(c * a[1][0]));
      REQUIRE((long long)(c * r[1][1]) == (long long)(c * a[1][1]));
      REQUIRE((long long)(c * r[1][2]) == (long long)(c * a[1][2]));
      REQUIRE((long long)(c * r[2][0]) == (long long)(c * a[2][0]));
      REQUIRE((long long)(c * r[2][1]) == (long long)(c * a[2][1]));
      REQUIRE((long long)(c * r[2][2]) == (long long)(c * a[2][2]));
   }

   SECTION("triclinic2")
   {
      double m[3][3] = {{1, 1. / 3., 0.5}, //
                        {0, 2, 0.2},       //
                        {0, 0, 1}};
      double r[3][3] = {{1.105170918, 0.03874394669, 0.05563952912}, //
                        {0, 1.221402758, 0.02324636802},             //
                        {0, 0, 1.105170918}};

      double a[3][3] = {0};
      trimat_exp(a, m, 0.1);

      REQUIRE((long long)(c * r[0][0]) == (long long)(c * a[0][0]));
      REQUIRE((long long)(c * r[0][1]) == (long long)(c * a[0][1]));
      REQUIRE((long long)(c * r[0][2]) == (long long)(c * a[0][2]));
      REQUIRE((long long)(c * r[1][0]) == (long long)(c * a[1][0]));
      REQUIRE((long long)(c * r[1][1]) == (long long)(c * a[1][1]));
      REQUIRE((long long)(c * r[1][2]) == (long long)(c * a[1][2]));
      REQUIRE((long long)(c * r[2][0]) == (long long)(c * a[2][0]));
      REQUIRE((long long)(c * r[2][1]) == (long long)(c * a[2][1]));
      REQUIRE((long long)(c * r[2][2]) == (long long)(c * a[2][2]));
   }
}
