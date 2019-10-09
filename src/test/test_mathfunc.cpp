#include "mathfunc.h"
#include "test.h"

using namespace TINKER_NAMESPACE;

TEST_CASE("MathFunc", "[util]")
{
   SECTION("IsPow2")
   {
      auto f = is_pow2;

      REQUIRE(f(0) == false);
      REQUIRE(f(1) == true);
      REQUIRE(f(2) == true);
      REQUIRE(f(3) == false);
      REQUIRE(f(4) == true);
      REQUIRE(f(5) == false);
      REQUIRE(f(6) == false);
      REQUIRE(f(7) == false);
      REQUIRE(f(8) == true);
   }

   SECTION("IntLog2")
   {
      auto f = int_log2;

      REQUIRE(f(1) == 0);
      REQUIRE(f(2) == 1);
      REQUIRE(f(3) == 1);
      REQUIRE(f(4) == 2);
      REQUIRE(f(5) == 2);
      REQUIRE(f(6) == 2);
      REQUIRE(f(7) == 2);
      REQUIRE(f(8) == 3);
   }

   SECTION("Pow2LessOrEqual")
   {
      auto f = pow2_le;

      REQUIRE(f(1) == 1);
      REQUIRE(f(2) == 2);
      REQUIRE(f(3) == 2);
      REQUIRE(f(4) == 4);
      REQUIRE(f(5) == 4);
      REQUIRE(f(6) == 4);
      REQUIRE(f(7) == 4);
      REQUIRE(f(8) == 8);
      REQUIRE(f(9) == 8);
      REQUIRE(f(1023) == 512);
      REQUIRE(f(1024) == 1024);
      REQUIRE(f(1025) == 1024);
   }

   SECTION("Pow2GreaterOrEqual")
   {
      auto f = pow2_ge;

      REQUIRE(f(0) == 1);
      REQUIRE(f(1) == 1);
      REQUIRE(f(2) == 2);
      REQUIRE(f(3) == 4);
      REQUIRE(f(4) == 4);
      REQUIRE(f(5) == 8);
      REQUIRE(f(6) == 8);
      REQUIRE(f(7) == 8);
      REQUIRE(f(8) == 8);
      REQUIRE(f(9) == 16);
      REQUIRE(f(1023) == 1024);
      REQUIRE(f(1024) == 1024);
   }
}
