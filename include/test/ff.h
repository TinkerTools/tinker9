#ifndef TINKER_TEST_FF_H_
#define TINKER_TEST_FF_H_

#include "gpu/decl_mdstate.h"
#include "gpu/e_potential.h"
#include "gpu/gpu.h"
#include "util_text.h"
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
    REQUIRE(eng == Approx(ref_eng).margin(eps));                               \
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
        REQUIRE(vir[k] == Approx(ref_v[i][j]).margin(eps));                    \
      }                                                                        \
    }                                                                          \
  }
#define COMPARE_VIR2_(gpuptr, gpuptr2, ref_v, eps)                             \
  {                                                                            \
    double vir[9], vir2[9];                                                    \
    gpu::get_virial(vir, gpuptr);                                              \
    gpu::get_virial(vir2, gpuptr2);                                            \
    for (int i = 0; i < 3; ++i) {                                              \
      for (int j = 0; j < 3; ++j) {                                            \
        int k = 3 * i + j;                                                     \
        REQUIRE((vir[k] + vir2[k]) == Approx(ref_v[i][j]).margin(eps));        \
      }                                                                        \
    }                                                                          \
  }
#define COMPARE_GRADIENT3_(gx, gy, gz, ref_grad, eps, check_ij)                \
  {                                                                            \
    grad_t grad(gpu::n);                                                       \
    double* dst = &grad[0][0];                                                 \
    gpu::copyout_array2(0, 3, dst, gx, gpu::n);                                \
    gpu::copyout_array2(1, 3, dst, gy, gpu::n);                                \
    gpu::copyout_array2(2, 3, dst, gz, gpu::n);                                \
    for (int i = 0; i < gpu::n; ++i) {                                         \
      for (int j = 0; j < 3; ++j) {                                            \
        if (check_ij(i, j))                                                    \
          REQUIRE(grad[i][j] == Approx(ref_grad[i][j]).margin(eps));           \
      }                                                                        \
    }                                                                          \
  }
#define COMPARE_GRADIENT2_(ref_grad, eps, check_ij)                            \
  COMPARE_GRADIENT3_(gpu::gx, gpu::gy, gpu::gz, ref_grad, eps, check_ij)
#define COMPARE_GRADIENT_(ref_grad, eps)                                       \
  COMPARE_GRADIENT2_(ref_grad, eps, [](int, int) { return true; })
#define PRINT_ENERGY_(gpuptr)                                                  \
  {                                                                            \
    double e_ = gpu::get_energy(gpuptr);                                       \
    print(stdout, " ENERGY{:>12.4f}\n", e_);                                   \
  }
#define PRINT_COUNT_(gpuptr)                                                   \
  {                                                                            \
    int c_ = gpu::get_count(gpuptr);                                           \
    print(stdout, " COUNT{:>12d}\n", c_);                                      \
  }
#define PRINT_VIR_(gpuptr)                                                     \
  {                                                                            \
    double vir[9];                                                             \
    gpu::get_virial(vir, gpuptr);                                              \
    for (int i = 0; i < 3; ++i) {                                              \
      print(stdout, " VIRIAL{:>12.4f}{:>12.4f}{:>12.4f}\n", vir[3 * i],        \
            vir[3 * i + 1], vir[3 * i + 2]);                                   \
    }                                                                          \
  }
#define PRINT_VIR2_(gpuptr, gpuptr2)                                           \
  {                                                                            \
    double vir[9], vir2[9];                                                    \
    gpu::get_virial(vir, gpuptr);                                              \
    gpu::get_virial(vir2, gpuptr2);                                            \
    for (int i = 0; i < 9; ++i)                                                \
      vir[i] += vir2[i];                                                       \
    for (int i = 0; i < 3; ++i) {                                              \
      print(stdout, " VIRIAL{:>12.4f}{:>12.4f}{:>12.4f}\n", vir[3 * i],        \
            vir[3 * i + 1], vir[3 * i + 2]);                                   \
    }                                                                          \
  }
#define PRINT_V3_(info, v1, v2, v3)                                            \
  {                                                                            \
    grad_t grad(gpu::n);                                                       \
    double* dst = &grad[0][0];                                                 \
    gpu::copyout_array2(0, 3, dst, v1, gpu::n);                                \
    gpu::copyout_array2(1, 3, dst, v2, gpu::n);                                \
    gpu::copyout_array2(2, 3, dst, v3, gpu::n);                                \
    for (int i = 0; i < gpu::n; ++i) {                                         \
      print(stdout, " {:s} ATOM{:>6d}{:>12.6f}{:>12.6f}{:>12.6f}\n", info,     \
            i + 1, grad[i][0], grad[i][1], grad[i][2]);                        \
    }                                                                          \
  }
#define PRINT_GRADIENT_ PRINT_V3_("GRADIENT", gpu::gx, gpu::gy, gpu::gz)
#define PRINT_COORD_ PRINT_V3_("COORD", gpu::x, gpu::y, gpu::z)
#define PRINT_VELOCITY_ PRINT_V3_("VELOCITY", gpu::vx, gpu::vy, gpu::vz)
#define PRINT_GRAD_T_(grad)                                                    \
  {                                                                            \
    for (int i = 0; i < grad.size(); ++i) {                                    \
      print(stdout, " ATOM{:>6d}{:>12.4f}{:>12.4f}{:>12.4f}\n", i + 1,         \
            grad[i][0], grad[i][1], grad[i][2]);                               \
    }                                                                          \
  }
#define COMPARE_BONED_FORCE(routine, gpu_e, ref_e, eps_e, cpu_count,           \
                            ref_count, gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g,   \
                            gpu_v, ref_v, eps_v)                               \
  {                                                                            \
    auto do_ij_ = [](int, int) { return true; };                               \
    gpu::zero_egv();                                                           \
    routine(gpu::v3);                                                          \
    COMPARE_ENERGY_(gpu_e, ref_e, eps_e);                                      \
    REQUIRE(cpu_count == ref_count);                                           \
                                                                               \
    gpu::zero_egv();                                                           \
    routine(gpu::v1);                                                          \
    COMPARE_ENERGY_(gpu_e, ref_e, eps_e);                                      \
    COMPARE_GRADIENT3_(gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g, do_ij_);          \
    COMPARE_VIR_(gpu_v, ref_v, eps_v);                                         \
                                                                               \
    gpu::zero_egv();                                                           \
    routine(gpu::v4);                                                          \
    COMPARE_ENERGY_(gpu_e, ref_e, eps_e);                                      \
    COMPARE_GRADIENT3_(gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g, do_ij_);          \
                                                                               \
    gpu::zero_egv();                                                           \
    routine(gpu::v5);                                                          \
    COMPARE_GRADIENT3_(gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g, do_ij_);          \
                                                                               \
    gpu::zero_egv();                                                           \
    routine(gpu::v6);                                                          \
    COMPARE_GRADIENT3_(gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g, do_ij_);          \
    COMPARE_VIR_(gpu_v, ref_v, eps_v);                                         \
  }
TINKER_NAMESPACE_END

#endif
