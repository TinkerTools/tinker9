#include "test.h"
#include "tool/orthomatrix.h"
using namespace tinker;

/* //
static auto prteig = [](double w[3]) {
   printf("%12.4lf%12.4lf%12.4lf\n", w[0], w[1], w[2]);
};
static auto prtmat = [](double m[3][3]) {
   printf("%12.4lf%12.4lf%12.4lf\n", m[0][0], m[0][1], m[0][2]);
   printf("%12.4lf%12.4lf%12.4lf\n", m[1][0], m[1][1], m[1][2]);
   printf("%12.4lf%12.4lf%12.4lf\n", m[2][0], m[2][1], m[2][2]);
};
// */

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
}
