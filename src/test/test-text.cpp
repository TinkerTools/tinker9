#include "test/test.h"
#include "util/text.h"

m_tinker_using_namespace;

TEST_CASE("Text") {
  SECTION("Replace") {
    char c = ' ';
    std::string s, r, ans;
    auto f = Text::replace;

    s = "";
    r = "x";
    ans = "";
    f(s, r, c);
    REQUIRE(s == ans);

    s = "1234567890";
    r = "42";
    ans = "1 3 567890";
    f(s, r, c);
    REQUIRE(s == ans);
  }
}
