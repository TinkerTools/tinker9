#ifndef TINKER_UTIL_TEST_RT_H_
#define TINKER_UTIL_TEST_RT_H_

#include "array.h"
#include "md.h"
#include "rc_man.h"
#include "util_io.h"
#include "util_potential.h"
#include <array>
#include <string>
#include <vector>

TINKER_NAMESPACE_BEGIN
/**
 * @brief
 * write a file to the disk in its constructor and remove this file in its
 * destructor
 */
class TestFile {
private:
  bool good_;
  std::string name_;

public:
  TestFile(const std::string& name, const std::string& content);
  ~TestFile();
};

/**
 * @brief
 * if possible, remove the file with the given name in its destructor
 */
class TestFileExpected {
private:
  std::string name_;

public:
  TestFileExpected(const std::string& name);
  ~TestFileExpected();
};

double test_get_eps(double eps_single, double eps_double);

void test_begin_with_args(int argc, const char** argv);
void test_end();
void test_mdinit(double t = 0, double atm = 0);
TINKER_NAMESPACE_END

#define COMPARE_ENERGY_(gpuptr, ref_eng, eps)                                  \
  {                                                                            \
    double eng = get_energy(gpuptr);                                           \
    REQUIRE(eng == Approx(ref_eng).margin(eps));                               \
  }
#define COMPARE_COUNT_(gpuptr, ref_count)                                      \
  {                                                                            \
    int count = get_count(gpuptr);                                             \
    REQUIRE(count == ref_count);                                               \
  }
#define COMPARE_VIR_(gpuptr, ref_v, eps)                                       \
  {                                                                            \
    double vir1[9];                                                            \
    get_virial(vir1, gpuptr);                                                  \
    for (int i = 0; i < 3; ++i) {                                              \
      for (int j = 0; j < 3; ++j) {                                            \
        int k = 3 * i + j;                                                     \
        REQUIRE(vir1[k] == Approx(ref_v[i][j]).margin(eps));                   \
      }                                                                        \
    }                                                                          \
  }
#define COMPARE_VIR2_(gpuptr, gpuptr2, ref_v, eps)                             \
  {                                                                            \
    double vir1[9], vir2[9];                                                   \
    get_virial(vir1, gpuptr);                                                  \
    get_virial(vir2, gpuptr2);                                                 \
    for (int i = 0; i < 3; ++i) {                                              \
      for (int j = 0; j < 3; ++j) {                                            \
        int k = 3 * i + j;                                                     \
        REQUIRE((vir1[k] + vir2[k]) == Approx(ref_v[i][j]).margin(eps));       \
      }                                                                        \
    }                                                                          \
  }
#define COMPARE_GRADIENT3_(gx, gy, gz, ref_grad, eps, check_ij)                \
  {                                                                            \
    std::vector<std::array<double, 3>> grad(n);                                \
    double* dst = &grad[0][0];                                                 \
    copyout_array2(0, 3, dst, gx, n);                                          \
    copyout_array2(1, 3, dst, gy, n);                                          \
    copyout_array2(2, 3, dst, gz, n);                                          \
    for (int i = 0; i < n; ++i) {                                              \
      for (int j = 0; j < 3; ++j) {                                            \
        if (check_ij(i, j))                                                    \
          REQUIRE(grad[i][j] == Approx(ref_grad[i][j]).margin(eps));           \
      }                                                                        \
    }                                                                          \
  }
#define COMPARE_GRADIENT2_(ref_grad, eps, check_ij)                            \
  COMPARE_GRADIENT3_(gx, gy, gz, ref_grad, eps, check_ij)
#define COMPARE_GRADIENT_(ref_grad, eps)                                       \
  COMPARE_GRADIENT2_(ref_grad, eps, [](int, int) { return true; })
#define COMPARE_BONDED_FORCE(routine, gpu_e, ref_e, eps_e, cpu_count,          \
                             ref_count, gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g,  \
                             gpu_v, ref_v, eps_v)                              \
  {                                                                            \
    auto do_ij_ = [](int, int) { return true; };                               \
    zero_egv();                                                                \
    routine(calc::v3);                                                         \
    COMPARE_ENERGY_(gpu_e, ref_e, eps_e);                                      \
    REQUIRE(cpu_count == ref_count);                                           \
                                                                               \
    zero_egv();                                                                \
    routine(calc::v1);                                                         \
    COMPARE_ENERGY_(gpu_e, ref_e, eps_e);                                      \
    COMPARE_GRADIENT3_(gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g, do_ij_);          \
    COMPARE_VIR_(gpu_v, ref_v, eps_v);                                         \
                                                                               \
    zero_egv();                                                                \
    routine(calc::v4);                                                         \
    COMPARE_ENERGY_(gpu_e, ref_e, eps_e);                                      \
    COMPARE_GRADIENT3_(gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g, do_ij_);          \
                                                                               \
    zero_egv();                                                                \
    routine(calc::v5);                                                         \
    COMPARE_GRADIENT3_(gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g, do_ij_);          \
                                                                               \
    zero_egv();                                                                \
    routine(calc::v6);                                                         \
    COMPARE_GRADIENT3_(gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g, do_ij_);          \
    COMPARE_VIR_(gpu_v, ref_v, eps_v);                                         \
  }

#endif
