#include "box.h"
#include "mathfunc.h"
#include "seq_image.h"
#include "test.h"
#include "tinker_rt.h"
#include <cmath>
#include <tinker/detail/boxes.hh>
using namespace TINKER_NAMESPACE;


static void set_box(Box& b, const double* p)
{
   boxes::orthogonal = 0;
   boxes::monoclinic = 0;
   boxes::triclinic = 0;
   boxes::octahedron = 0;
   if (b.shape == Box::ortho)
      boxes::orthogonal = 1;
   else if (b.shape == Box::mono)
      boxes::monoclinic = 1;
   else if (b.shape == Box::tri)
      boxes::triclinic = 1;

   boxes::xbox = p[0];
   boxes::ybox = p[1];
   boxes::zbox = p[2];
   boxes::alpha = p[3];
   boxes::beta = p[4];
   boxes::gamma = p[5];
   TINKER_RT(lattice)();

   b.volbox = boxes::volbox;
   real t[3][3], u[3][3];
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         t[i][j] = boxes::lvec[i][j];
         u[i][j] = boxes::recip[i][j];
      }
   }
   std::memcpy(&b.lvec[0][0], &t[0][0], sizeof(real) * 9);
   std::memcpy(&b.recip[0][0], &u[0][0], sizeof(real) * 9);
}


#define compare_im()                                                           \
   {                                                                           \
      real xx = xr, yy = yr, zz = zr;                                          \
      image_general(xx, yy, zz, &b);                                           \
      REQUIRE(xx == Approx(xa).margin(eps));                                   \
      REQUIRE(yy == Approx(ya).margin(eps));                                   \
      REQUIRE(zz == Approx(za).margin(eps));                                   \
   }
#define compare_in()                                                           \
   {                                                                           \
      real xx = xr, yy = yr, zz = zr;                                          \
      imagen_general(xx, yy, zz, &b);                                          \
      REQUIRE(REAL_ABS(xx) == Approx(std::fabs(xa)).margin(eps));              \
      REQUIRE(REAL_ABS(yy) == Approx(std::fabs(ya)).margin(eps));              \
      REQUIRE(REAL_ABS(zz) == Approx(std::fabs(za)).margin(eps));              \
   }
#define COS(x) std::cos(x* _1radian)
#define SIN(x) std::sin(x* _1radian)

TEST_CASE("Box-1", "[ff][box][orthogonal]")
{
   const char* argv[] = {"dummy"};
   int argc = 1;
   Box b;
   real xr, yr, zr;

   double eps = 1.0e-6;
   b.shape = Box::ortho;
   double p[] = {16, 16, 16, 90, 90, 90};
   set_box(b, p);

   fortran_runtime_initialize(argc, (char**)argv);
   TINKER_RT(initial)();

   SECTION("  - volume")
   {
      double vol = p[0] * p[1] * p[2];
      REQUIRE(b.volbox == Approx(vol).margin(eps));
   }

   // image and imagen
   SECTION("  - -1 to -1/2")
   {
      xr = -16, yr = -12, zr = -8;
      double xa = 0, ya = 4, za = -8;
      compare_im();
      compare_in();
   }

   SECTION("  - -1/4 to 1/4")
   {
      xr = -4, yr = 0, zr = 4;
      double xa = -4, ya = 0, za = 4;
      compare_im();
      compare_in();
   }

   SECTION("  - 1/2 to 1")
   {
      xr = 8, yr = 12, zr = 16;
      double xa = -8, ya = -4, za = 0;
      compare_im();
      compare_in();
   }

   TINKER_RT(final)();
   fortran_runtime_finish();
}


TEST_CASE("Box-2", "[ff][box][monoclinic]")
{
   const char* argv[] = {"dummy"};
   int argc = 1;
   Box b;
   real xr, yr, zr;

   double eps = 1.0e-6;
   b.shape = Box::mono;
   double p[] = {32, 24, 20, 90, 30, 90};
   set_box(b, p);

   fortran_runtime_initialize(argc, (char**)argv);
   TINKER_RT(initial)();

   SECTION(" -- lvec")
   {
      REQUIRE(b.lvec[0][0] == Approx(p[0]).margin(eps));
      REQUIRE(b.lvec[1][0] == Approx(0).margin(eps));
      REQUIRE(b.lvec[2][0] == Approx(0).margin(eps));

      REQUIRE(b.lvec[0][1] == Approx(0).margin(eps));
      REQUIRE(b.lvec[1][1] == Approx(p[1]).margin(eps));
      REQUIRE(b.lvec[2][1] == Approx(0).margin(eps));

      REQUIRE(b.lvec[0][2] == Approx(p[2] * COS(p[4])).margin(eps));
      REQUIRE(b.lvec[1][2] == Approx(0).margin(eps));
      REQUIRE(b.lvec[2][2] == Approx(p[2] * SIN(p[4])).margin(eps));
   }

   SECTION("  - volume")
   {
      double vol = p[0] * p[1] * p[2] * SIN(p[4]);
      REQUIRE(b.volbox == Approx(vol).margin(eps));
   }

   // image and imagen
   SECTION("  - origin")
   {
      xr = 0, yr = 0, zr = 0;
      double xa = 0, ya = 0, za = 0;
      compare_im();
      compare_in();
   }

   SECTION(" -- xy-plane")
   {
      xr = -8, yr = -6, zr = 0;
      double xa = -8, ya = -6, za = 0;
      compare_im();
      compare_in();
   }

   SECTION(" -- general")
   {
      double xa, ya, za;

      xr = 5, yr = 10, zr = 15;            // OA
      xa = 2.3589838486, ya = 10, za = -5; // OA'
      compare_im();
      compare_in();

      xr = -13, yr = -30, zr = 20;            // OB
      xa = -15.6410161514, ya = -6.0, za = 0; // OB'
      compare_im();
      compare_in();

      xr = -18, yr = -40, zr = 5; // AB = OB - OA
      xa = -3.3205080757, ya = 8.0, za = -5;
      compare_im();
      compare_in();

      xr = -18, yr = -16, zr = 5; // AB' = OB' - OA'
      xa = -3.3205080757, ya = 8.0, za = -5;
      compare_im();
      compare_in();
   }

   TINKER_RT(final)();
   fortran_runtime_finish();
}


TEST_CASE("Box-3", "[ff][box][triclinic]")
{
   const char* argv[] = {"dummy"};
   int argc = 1;
   Box b;
   real xr, yr, zr;
   double xa, ya, za;

   double eps = 1.0e-6;
   b.shape = Box::tri;
   double p[] = {32, 24, 20, 75, 60, 45};
   set_box(b, p);

   fortran_runtime_initialize(argc, (char**)argv);
   TINKER_RT(initial)();

   // volume
   double al = COS(p[3]);
   double be = COS(p[4]);
   double ga = COS(p[5]);
   double sq = 1.0 - al * al - be * be - ga * ga + 2 * al * be * ga;
   double vol = p[0] * p[1] * p[2] * std::sqrt(sq);
   REQUIRE(b.volbox == Approx(vol).margin(eps));

   // image and imagen
   SECTION(" -- origin")
   {
      xr = 0, yr = 0, zr = 0;
      xa = 0, ya = 0, za = 0;
      compare_im();
      compare_in();
   }

   SECTION(" -- general")
   {
      double xa, ya, za;

      xr = 5, yr = 10, zr = 15;                             // OA
      xa = 10.02943725, ya = -4.29107082, za = -2.11199354; // OA'
      compare_im();
      compare_in();

      xr = -13, yr = -30, zr = 20;                        // OB
      xa = 10.94112550, ya = 6.62061742, za = 2.88800646; // OB'
      compare_im();
      compare_in();

      xr = -18, yr = -40, zr = 5; // AB = OB - OA
      xa = -16.05887450, ya = -6.05887450, za = 5;
      compare_im();
      compare_in();

      xr = 0.91168825, yr = 10.91168824, zr = 5; // AB' = OB' - OA'
      xa = -16.05887450, ya = -6.05887450, za = 5;
      compare_im();
      compare_in();
   }

   TINKER_RT(final)();
   fortran_runtime_finish();
}
