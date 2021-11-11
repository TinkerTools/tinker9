#include "test.h"
#include "test_rt.h"
#include "tool/fc.h"


TEST_CASE("Tinker_RT", "[util]")
{
   SECTION("Version")
   {
      std::string ans, ref;

      ans = tinker_f_version(TINKER9_DIRSTR "/CMakeLists.txt", "old");
      ref = TINKER9_DIRSTR "/CMakeLists.txt";
      REQUIRE(ans == ref);

      ans = tinker_f_version(TINKER9_DIRSTR "/CMakeLists.txt", "new");
      ref = TINKER9_DIRSTR "/CMakeLists.txt_2";
      REQUIRE(ans == ref);
   }
}
