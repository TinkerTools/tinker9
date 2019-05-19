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

#define COMPARE_ENERGY_(gpuptr, ref_eng, eps)                                  \
  {                                                                            \
    double eng = gpu::get_energy(gpuptr);                                      \
    REQUIRE(eng == Approx(ref_eng).epsilon(eps));                              \
  }
#define COMPARE_COUNT_(gpuptr, ref_count)                                      \
  {                                                                            \
    int count = gpu::get_count(gpuptr);                                        \
    REQUIRE(count == ref_count);                                               \
  }
#define COMPARE_VIR_(gpuptr, ref_v, eps)                                       \
  {                                                                            \
    double vir[9];                                                             \
    gpu::get_virial(vir, gpuptr);                                              \
    for (int i = 0; i < 3; ++i) {                                              \
      for (int j = 0; j < 3; ++j) {                                            \
        int k = 3 * i + j;                                                     \
        REQUIRE(vir[k] == Approx(ref_v[i][j]).epsilon(eps));                   \
      }                                                                        \
    }                                                                          \
  }

TINKER_NAMESPACE_END

#endif
