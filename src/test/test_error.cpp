#define TINKER_ALWAYS_CHECK_CUDART

#include "test/test.h"
#include "util/error.h"
#include "util/error_cudart.h"
#include "util/format_print.h"

m_tinker_using_namespace;

TEST_CASE("Error", "[noassert][util]") {
  SECTION("TinkerThrowMacro") {
    try {
      m_tinker_throw("TinkerThrowMacro Test");
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }

  SECTION("CudaRTError1") {
    auto func = []() { return cudaErrorInvalidHostPointer; };
    try {
      check_cudart(func());
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }

  SECTION("CudaRTError2") {
    auto func = []() { return cudaErrorInvalidHostPointer; };
    try {
      check_cudart(func(), "Optional Error Message.");
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }

  SECTION("CudaRTError3") {
    auto func = []() { return cudaErrorInvalidHostPointer; };
    try {
      check_cudart(func(), cudaError_t, cudaSuccess);
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }
}
