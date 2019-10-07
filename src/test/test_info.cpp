#include "macro.h"
#include "test.h"

TINKER_NAMESPACE_BEGIN
extern void x_info (int, char**);
TINKER_NAMESPACE_END

using namespace TINKER_NAMESPACE;

TEST_CASE ("Info", "[noassert][info]")
{
   x_info (0, nullptr);
}
