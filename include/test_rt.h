#pragma once
#include "energy.h"
#include "md.h"
#include "rc_man.h"
#include <array>
#include <string>
#include <vector>


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup test
 * \brief Write a file to the disk in its constructor and remove this file in
 * its destructor, unless the file is set to be kept on disk.
 */
class TestFile
{
private:
   bool good_;
   std::string name_;

public:
   TestFile(const std::string& name, const std::string& content);
   ~TestFile();
   void keep();
};


/**
 * \ingroup
 * \brief Remove the file with the given name in its destructor if possible.
 */
class TestFileExpected
{
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
   {                                                                           \
      double eng = get_energy(gpuptr);                                         \
      REQUIRE(eng == Approx(ref_eng).margin(eps));                             \
   }
#define COMPARE_COUNT_(gpuptr, ref_count)                                      \
   {                                                                           \
      int count = get_count(gpuptr);                                           \
      REQUIRE(count == ref_count);                                             \
   }
#define COMPARE_VIR_(gpuptr, ref_v, eps)                                       \
   {                                                                           \
      real vir1[9];                                                            \
      get_virial(vir1, gpuptr);                                                \
      REQUIRE(vir1[0] == Approx(ref_v[0][0]).margin(eps));                     \
      REQUIRE(vir1[1] == Approx(ref_v[0][1]).margin(eps));                     \
      REQUIRE(vir1[2] == Approx(ref_v[0][2]).margin(eps));                     \
      REQUIRE(vir1[3] == Approx(ref_v[1][0]).margin(eps));                     \
      REQUIRE(vir1[4] == Approx(ref_v[1][1]).margin(eps));                     \
      REQUIRE(vir1[5] == Approx(ref_v[1][2]).margin(eps));                     \
      REQUIRE(vir1[6] == Approx(ref_v[2][0]).margin(eps));                     \
      REQUIRE(vir1[7] == Approx(ref_v[2][1]).margin(eps));                     \
      REQUIRE(vir1[8] == Approx(ref_v[2][2]).margin(eps));                     \
   }
#define COMPARE_VIR2_(gpuptr, gpuptr2, ref_v, eps)                             \
   {                                                                           \
      real vir1[9], vir2[9];                                                   \
      get_virial(vir1, gpuptr);                                                \
      get_virial(vir2, gpuptr2);                                               \
      REQUIRE(vir1[0] + vir2[0] == Approx(ref_v[0][0]).margin(eps));           \
      REQUIRE(vir1[1] + vir2[1] == Approx(ref_v[0][1]).margin(eps));           \
      REQUIRE(vir1[2] + vir2[2] == Approx(ref_v[0][2]).margin(eps));           \
      REQUIRE(vir1[3] + vir2[3] == Approx(ref_v[1][0]).margin(eps));           \
      REQUIRE(vir1[4] + vir2[4] == Approx(ref_v[1][1]).margin(eps));           \
      REQUIRE(vir1[5] + vir2[5] == Approx(ref_v[1][2]).margin(eps));           \
      REQUIRE(vir1[6] + vir2[6] == Approx(ref_v[2][0]).margin(eps));           \
      REQUIRE(vir1[7] + vir2[7] == Approx(ref_v[2][1]).margin(eps));           \
      REQUIRE(vir1[8] + vir2[8] == Approx(ref_v[2][2]).margin(eps));           \
   }
#define COMPARE_GRADIENT3_(gx, gy, gz, ref_grad, eps, check_ij)                \
   {                                                                           \
      std::vector<std::array<double, 3>> grad(n);                              \
      double* dst = &grad[0][0];                                               \
      device_array::copyout2(false, 0, 3, n, dst, gx);                         \
      device_array::copyout2(false, 1, 3, n, dst, gy);                         \
      device_array::copyout2(false, 2, 3, n, dst, gz);                         \
      for (int i = 0; i < n; ++i) {                                            \
         for (int j = 0; j < 3; ++j) {                                         \
            if (check_ij(i, j))                                                \
               REQUIRE(grad[i][j] == Approx(ref_grad[i][j]).margin(eps));      \
         }                                                                     \
      }                                                                        \
   }
#define COMPARE_GRADIENT2_(ref_grad, eps, check_ij)                            \
   COMPARE_GRADIENT3_(gx, gy, gz, ref_grad, eps, check_ij)
#define COMPARE_GRADIENT_(ref_grad, eps)                                       \
   COMPARE_GRADIENT2_(ref_grad, eps, [](int, int) { return true; })
#define COMPARE_BONDED_FORCE(routine, cpu_count, gpu_e, gpu_v, ref_e, eps_e,   \
                             ref_count, gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g,  \
                             ref_v, eps_v)                                     \
   {                                                                           \
      auto do_ij_ = [](int, int) { return true; };                             \
      zero_egv();                                                              \
      routine(calc::v3);                                                       \
      COMPARE_ENERGY_(gpu_e, ref_e, eps_e);                                    \
      REQUIRE(cpu_count == ref_count);                                         \
                                                                               \
      zero_egv();                                                              \
      routine(calc::v1);                                                       \
      COMPARE_ENERGY_(gpu_e, ref_e, eps_e);                                    \
      COMPARE_GRADIENT3_(gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g, do_ij_);        \
      COMPARE_VIR_(gpu_v, ref_v, eps_v);                                       \
                                                                               \
      zero_egv();                                                              \
      routine(calc::v4);                                                       \
      COMPARE_ENERGY_(gpu_e, ref_e, eps_e);                                    \
      COMPARE_GRADIENT3_(gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g, do_ij_);        \
                                                                               \
      zero_egv();                                                              \
      routine(calc::v5);                                                       \
      COMPARE_GRADIENT3_(gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g, do_ij_);        \
                                                                               \
      zero_egv();                                                              \
      routine(calc::v6);                                                       \
      COMPARE_GRADIENT3_(gpu_gx, gpu_gy, gpu_gz, ref_g, eps_g, do_ij_);        \
      COMPARE_VIR_(gpu_v, ref_v, eps_v);                                       \
   }
