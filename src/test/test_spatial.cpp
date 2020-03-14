#include "box.h"
#include "seq_spatial_box.h"
#include "spatial.h"
#include "test.h"
using namespace TINKER_NAMESPACE;


template <int V>
static int test_frac_to_box(int px, int py, int pz, real fx, real fy, real fz)
{
   int ix, iy, iz;
   frac_to_ixyz(ix, iy, iz, px, py, pz, fx, fy, fz);
   if (V == 2)
      return spatial_v2::ixyz_to_box(px, py, pz, ix, iy, iz);
   else
      return spatial_v1::ixyz_to_box(px, py, pz, ix, iy, iz);
}


TEST_CASE("Spatial-V1", "[ff][spatial]")
{
   SECTION("  - 432")
   {
      int a[16][8][4];
      int* a0 = &a[0][0][0];
      int* a1 = &a[11][5][3];
      real fx = 11.5f / 16 - 0.5f;
      real fy = 5.5f / 8 - 0.5f;
      real fz = 3.5f / 4 - 0.5f;
      int id = test_frac_to_box<1>(4, 3, 2, fx, fy, fz);
      REQUIRE(id == (a1 - a0));


      int ix, iy, iz;
      spatial_v1::box_to_ixyz(ix, iy, iz, 4, 3, 2, id);
      REQUIRE(ix == 11);
      REQUIRE(iy == 5);
      REQUIRE(iz == 3);
   }
}


TEST_CASE("Spatial-V2", "[ff][spatial]")
{
   SECTION("  - 333")
   {
      // [8][8][8]
      // [5][3][2] [101][011][010] [100 011 110] [0x11E]
      real fx = 5.5f / 8 - 0.5f;
      real fy = 3.3f / 8 - 0.5f;
      real fz = 2.2f / 8 - 0.5f;
      int id = test_frac_to_box<2>(3, 3, 3, fx, fy, fz);
      REQUIRE(id == 0x11E);

      int ix, iy, iz;
      spatial_v2::box_to_ixyz(ix, iy, iz, 3, 3, 3, id);
      REQUIRE(ix == 5);
      REQUIRE(iy == 3);
      REQUIRE(iz == 2);
   }


   SECTION("  - 433")
   {
      // [16][8][8]
      // [13][3][2] [1101][011][010] [1 100 011 110] [0x31E]
      real fx = 13.6f / 16 - 0.5f;
      real fy = 3.3f / 8 - 0.5f;
      real fz = 2.2f / 8 - 0.5f;
      int id = test_frac_to_box<2>(4, 3, 3, fx, fy, fz);
      REQUIRE(id == 0x31E);

      int ix, iy, iz;
      spatial_v2::box_to_ixyz(ix, iy, iz, 4, 3, 3, id);
      REQUIRE(ix == 13);
      REQUIRE(iy == 3);
      REQUIRE(iz == 2);
   }


   SECTION("  - 443")
   {
      // [16][16][8]
      // [11][13][5] [1011][1101][101] [11 011 100 111] [0x6E7]
      real fx = 11.7f / 16 - 0.5f;
      real fy = 13.4f / 16 - 0.5f;
      real fz = 5.2f / 8 - 0.5f;
      int id = test_frac_to_box<2>(4, 4, 3, fx, fy, fz);
      REQUIRE(id == 0x6E7);

      int ix, iy, iz;
      spatial_v2::box_to_ixyz(ix, iy, iz, 4, 4, 3, id);
      REQUIRE(ix == 11);
      REQUIRE(iy == 13);
      REQUIRE(iz == 5);
   }
}


TEST_CASE("Spatial-Cut", "[ff][spatial]")
{
   lvec1 = make_real3(40, 0, .0);
   lvec2 = make_real3(0, 160, 0);
   lvec3 = make_real3(0, 0, 40);
   const int level = 4;
   int px, py, pz;


   SECTION("   - spatial_cut_v1")
   {
      spatial_cut_v1(px, py, pz, level);
      REQUIRE(px == 2);
      REQUIRE(py == 1);
      REQUIRE(pz == 1);
   }


   SECTION("   - spatial_cut_v2")
   {
      spatial_cut_v2(px, py, pz, level);
      REQUIRE(px == 1);
      REQUIRE(py == 1);
      REQUIRE(pz == 2);
   }


   SECTION("   - spatial_cut_v3")
   {
      spatial_cut_v3(px, py, pz, level);
      REQUIRE(px == 0);
      REQUIRE(py == 3);
      REQUIRE(pz == 1);
   }
}
