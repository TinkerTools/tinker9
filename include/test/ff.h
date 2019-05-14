#ifndef TINKER_TEST_FF_H_
#define TINKER_TEST_FF_H_

#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "gpu/e_potential.h"
#include <array>
#include <vector>

TINKER_NAMESPACE_BEGIN
namespace test {
/// @brief Gradient type used in unit tests.
typedef std::vector<std::array<double, 3>> grad_t;
}
TINKER_NAMESPACE_END

#endif
