#define TINKER_ALWAYS_CHECK_CUDART

#include "gpu/rc.h"
#include "test/test.h"
#include "util/format_print.h"

m_tinker_using_namespace;
using namespace gpu;

TEST_CASE("Error", "[noassert][util]") {
  const char* fmt = "=== end of this test section ===\n\n\n";

  SECTION("TinkerThrowMacro") {
    try {
      m_tinker_throw("TinkerThrowMacro Test");
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    print(stdout, fmt);
  }

  SECTION("CudaRTError1") {
    auto func = []() { return cudaErrorInvalidHostPointer; };
    try {
      check_cudart(func());
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    print(stdout, fmt);
  }

  SECTION("CudaRTError2") {
    auto func = []() { return cudaErrorInvalidHostPointer; };
    try {
      check_cudart(func(), "Optional Error Message.");
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    print(stdout, fmt);
  }

  SECTION("CudaRTError3") {
    auto func = []() { return cudaErrorInvalidHostPointer; };
    try {
      check_cudart(func(), cudaError_t, cudaSuccess);
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
    print(stdout, fmt);
  }
}
