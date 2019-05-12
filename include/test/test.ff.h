#ifndef TINKER_TEST_TEST_FF_H_
#define TINKER_TEST_TEST_FF_H_

#include "gpu/data.h"
#include "gpu/mdstate.h"
#include "gpu/potential.h"
#include <array>
#include <vector>

TINKER_NAMESPACE_BEGIN
namespace test {
typedef std::vector<std::array<double, 3>> grad_t;
}
TINKER_NAMESPACE_END

#endif
