#define TINKER_ALWAYS_CHECK_RT

#include "util_rt.h"
#include "util_test.h"

using namespace TINKER_NAMESPACE;

TEST_CASE("Error", "[noassert][error]") {
  const char* fmt = "=== end of this test section ===\n\n\n";

  SECTION("TinkerThrowMacro") {
    try {
      TINKER_THROW("TinkerThrowMacro Test");
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    print(stdout, fmt);
  }

  SECTION("CudaRTError1") {
    auto func = []() { return 42; };
    try {
      check_rt(func());
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    print(stdout, fmt);
  }

  SECTION("CudaRTError2") {
    auto func = []() { return 42; };
    try {
      check_rt(func(), "Optional Error Message");
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    print(stdout, fmt);
  }
}
