#include "test/test.h"
#include "util/error.h"
#include "util/format.print.h"

m_tinker_using_namespace;

TEST_CASE("UtilError", "[act]") {
  SECTION("TinkerThrowMacro") {
    try {
      m_tinker_throw("TinkerThrowMacro Test");
    } catch (std::exception& e) {
      std::cout << e.what() << std::endl;
    }
  }
}
