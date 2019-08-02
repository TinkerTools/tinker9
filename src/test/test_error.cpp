#define TINKER_ALWAYS_CHECK_CUDART

#include "test/test.h"
#include "util_rt.h"

using namespace TINKER_NAMESPACE;

TEST_CASE("Error", "[noassert][util]") {
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
    auto func = []() { return cudaErrorInvalidHostPointer; };
    try {
      check_rt(func());
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    print(stdout, fmt);
  }

  SECTION("CudaRTError2") {
    auto func = []() { return cudaErrorInvalidHostPointer; };
    try {
      check_rt(func(), "Optional Error Message.");
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    print(stdout, fmt);
  }

  SECTION("CudaRTError3") {
    auto func = []() { return cudaErrorInvalidHostPointer; };
    try {
      check_rt(func(), cudaError_t, cudaSuccess);
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    print(stdout, fmt);
  }
}
