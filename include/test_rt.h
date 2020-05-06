#pragma once
#include "energy.h"
#include "md.h"
#include "rc_man.h"
#include <array>
#include <string>
#include <vector>


namespace tinker {
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
}

#define COMPARE_INTS(i1, refi)       REQUIRE(i1 == refi)
#define COMPARE_REALS(v1, refv, eps) REQUIRE(v1 == Approx(refv).margin(eps))
#define COMPARE_ENERGY(gpuptr, ref_eng, eps)                                   \
   {                                                                           \
      double eng = energy_reduce(gpuptr);                                      \
      REQUIRE(eng == Approx(ref_eng).margin(eps));                             \
   }
#define COMPARE_COUNT(gpuptr, ref_count)                                       \
   {                                                                           \
      int count = count_reduce(gpuptr);                                        \
      REQUIRE(count == ref_count);                                             \
   }
#define COMPARE_VIR(gpuptr, ref_v, eps)                                        \
   {                                                                           \
      virial_prec vir1[9];                                                     \
      virial_reduce(vir1, gpuptr);                                             \
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
#define COMPARE_VIR2(gpuptr, gpuptr2, ref_v, eps)                              \
   {                                                                           \
      virial_prec vir1[9], vir2[9];                                            \
      virial_reduce(vir1, gpuptr);                                             \
      virial_reduce(vir2, gpuptr2);                                            \
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
#define COMPARE_GRADIENT2(ref_grad, eps, check_ij)                             \
   {                                                                           \
      std::vector<double> gradx(n), grady(n), gradz(n);                        \
      copy_gradient(calc::grad, gradx.data(), grady.data(), gradz.data());     \
      for (int i = 0; i < n; ++i) {                                            \
         if (check_ij(i, 0))                                                   \
            REQUIRE(gradx[i] == Approx(ref_grad[i][0]).margin(eps));           \
         if (check_ij(i, 1))                                                   \
            REQUIRE(grady[i] == Approx(ref_grad[i][1]).margin(eps));           \
         if (check_ij(i, 2))                                                   \
            REQUIRE(gradz[i] == Approx(ref_grad[i][2]).margin(eps));           \
      }                                                                        \
   }
#define COMPARE_GRADIENT(ref_grad, eps)                                        \
   COMPARE_GRADIENT2(ref_grad, eps, [](int, int) { return true; })
#define COMPARE_BONDED_FORCE(routine, cpu_count, gpu_e, gpu_v, ref_e, eps_e,   \
                             ref_count, ref_g, eps_g, ref_v, eps_v)            \
   {                                                                           \
      auto do_ij_ = [](int, int) { return true; };                             \
      zero_egv();                                                              \
      routine(calc::v3);                                                       \
      sum_energy(calc::v3);                                                    \
      COMPARE_ENERGY(gpu_e, ref_e, eps_e);                                     \
      REQUIRE(cpu_count == ref_count);                                         \
                                                                               \
      zero_egv();                                                              \
      routine(calc::v1);                                                       \
      sum_energy(calc::v1);                                                    \
      COMPARE_ENERGY(gpu_e, ref_e, eps_e);                                     \
      COMPARE_GRADIENT2(ref_g, eps_g, do_ij_);                                 \
      COMPARE_VIR(gpu_v, ref_v, eps_v);                                        \
                                                                               \
      zero_egv();                                                              \
      routine(calc::v4);                                                       \
      sum_energy(calc::v4);                                                    \
      COMPARE_ENERGY(gpu_e, ref_e, eps_e);                                     \
      COMPARE_GRADIENT2(ref_g, eps_g, do_ij_);                                 \
                                                                               \
      zero_egv();                                                              \
      routine(calc::v5);                                                       \
      sum_energy(calc::v5);                                                    \
      COMPARE_GRADIENT2(ref_g, eps_g, do_ij_);                                 \
                                                                               \
      zero_egv();                                                              \
      routine(calc::v6);                                                       \
      sum_energy(calc::v6);                                                    \
      COMPARE_GRADIENT2(ref_g, eps_g, do_ij_);                                 \
      COMPARE_VIR(gpu_v, ref_v, eps_v);                                        \
   }
