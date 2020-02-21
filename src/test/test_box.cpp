#include "box.h"
#include "mathfunc.h"
#include "seq_image.h"
#include "test.h"
#include "tinker_rt.h"
#include <cmath>
#include <tinker/detail/boxes.hh>
using namespace TINKER_NAMESPACE;


namespace {
void set_box(BoxShape shape, const double* p)
{
   boxes::orthogonal = 0;
   boxes::monoclinic = 0;
   boxes::triclinic = 0;
   boxes::octahedron = 0;
   if (shape == ORTHO_BOX)
      boxes::orthogonal = 1;
   else if (shape == MONO_BOX)
      boxes::monoclinic = 1;
   else if (shape == TRI_BOX)
      boxes::triclinic = 1;

   boxes::xbox = p[0];
   boxes::ybox = p[1];
   boxes::zbox = p[2];
   boxes::alpha = p[3];
   boxes::beta = p[4];
   boxes::gamma = p[5];
   TINKER_RT(lattice)();

   recipa =
      make_real3(boxes::recip[0][0], boxes::recip[0][1], boxes::recip[0][2]);
   recipb =
      make_real3(boxes::recip[1][0], boxes::recip[1][1], boxes::recip[1][2]);
   recipc =
      make_real3(boxes::recip[2][0], boxes::recip[2][1], boxes::recip[2][2]);
   lvec1 = make_real3(boxes::lvec[0][0], boxes::lvec[0][1], boxes::lvec[0][2]);
   lvec2 = make_real3(boxes::lvec[1][0], boxes::lvec[1][1], boxes::lvec[1][2]);
   lvec3 = make_real3(boxes::lvec[2][0], boxes::lvec[2][1], boxes::lvec[2][2]);
}
}


#define compare_im()                                                           \
   {                                                                           \
      real xx = xr, yy = yr, zz = zr;                                          \
      image(xx, yy, zz);                                                       \
      REQUIRE(xx == Approx(xa).margin(eps));                                   \
      REQUIRE(yy == Approx(ya).margin(eps));                                   \
      REQUIRE(zz == Approx(za).margin(eps));                                   \
   }
#define compare_in()                                                           \
   {                                                                           \
      real xx = xr, yy = yr, zz = zr;                                          \
      imagen2(xx, yy, zz);                                                     \
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
   real xr, yr, zr;

   double eps = 1.0e-6;
   double p[] = {16, 16, 16, 90, 90, 90};
   set_box(ORTHO_BOX, p);

   fortran_runtime_initialize(argc, (char**)argv);
   TINKER_RT(initial)();

   SECTION("  - volume")
   {
      double vol = p[0] * p[1] * p[2];
      REQUIRE(volbox() == Approx(vol).margin(eps));
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
   double p[] = {32, 24, 20, 90, 30, 90};
   set_box(MONO_BOX, p);

   fortran_runtime_initialize(argc, (char**)argv);
   TINKER_RT(initial)();

   SECTION(" -- lvec")
   {
      REQUIRE(lvec1.x == Approx(p[0]).margin(eps));
      REQUIRE(lvec2.x == Approx(0).margin(eps));
      REQUIRE(lvec3.x == Approx(0).margin(eps));

      REQUIRE(lvec1.y == Approx(0).margin(eps));
      REQUIRE(lvec2.y == Approx(p[1]).margin(eps));
      REQUIRE(lvec3.y == Approx(0).margin(eps));

      REQUIRE(lvec1.z == Approx(p[2] * COS(p[4])).margin(eps));
      REQUIRE(lvec2.z == Approx(0).margin(eps));
      REQUIRE(lvec3.z == Approx(p[2] * SIN(p[4])).margin(eps));
   }

   SECTION("  - volume")
   {
      double vol = p[0] * p[1] * p[2] * SIN(p[4]);
      REQUIRE(volbox() == Approx(vol).margin(eps));
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
   double p[] = {32, 24, 20, 75, 60, 45};
   set_box(TRI_BOX, p);

   fortran_runtime_initialize(argc, (char**)argv);
   TINKER_RT(initial)();

   // volume
   double al = COS(p[3]);
   double be = COS(p[4]);
   double ga = COS(p[5]);
   double sq = 1.0 - al * al - be * be - ga * ga + 2 * al * be * ga;
   double vol = p[0] * p[1] * p[2] * std::sqrt(sq);
   REQUIRE(volbox() == Approx(vol).margin(eps));

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
