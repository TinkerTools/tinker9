#pragma once
#include "energy.h"
#include "md.h"
#include "tool/rc_man.h"
#include <array>
#include <string>
#include <vector>


namespace tinker {
/**
 * \ingroup test
 * Writes a file to disk in its constructor and remove this file in
 * its destructor, unless the file is set to be kept.
 */
class TestFile
{
private:
   bool good;
   std::string name;


public:
   /**
    * Copies file from `src` to `dst` and append the `extra` text to `dst` if
    * `extra` is not empty.
    *
    * If `dst` is an empty string, `dst` will have the same filename in the
    * current working directory. If `src` is an empty string, it will create a
    * new file at `dst`.
    */
   TestFile(const std::string& src, std::string dst = "",
            std::string extra = "");
   /** Removes the file on disk if possible. */
   ~TestFile();
   /** Prevents the file being deleted. */
   void __keep();
};


/**
 * \ingroup test
 * Removes the file with the given name in its destructor if possible.
 */
class TestFileExpected
{
private:
   std::string name_;


public:
   /// \param name  Name of the file to delete.
   TestFileExpected(const std::string& name);
   ~TestFileExpected();
};


/**
 * \ingroup test
 * Read reference values from a text file.
 */
class TestReference
{
private:
   int count;
   double energy;
   double virial[3][3];
   std::vector<double> gradient;

public:
   TestReference(std::string txt_filename);
   int get_count() const;
   double get_energy() const;
   double (*get_virial())[3];
   double (*get_gradient())[3];
};


/**
 * \ingroup test
 * Returns tolerance eps depending on the predefined floating-point precision.
 * \param eps_single  Larger `eps` for lower floating-point precision.
 * \param eps_double  Smaller `eps` for higher floating-point precision.
 */
double test_get_eps(double eps_single, double eps_double);


/**
 * \ingroup test
 * Initializes the test.
 */
void test_begin_with_args(int argc, const char** argv);
/**
 * \ingroup test
 * Ends the test.
 */
void test_end();
/**
 * \ingroup test
 * Initializes MD in the test.
 * \param t    Temperature in Kelvin.
 * \param atm  Atmosphere in atm.
 */
void test_mdinit(double t = 0, double atm = 0);
}


/**
 * \def COMPARE_INTS
 * \ingroup test
 * Compare two integers.
 *
 * \def COMPARE_REALS
 * \ingroup test
 * Compare two floating-point numbers with a margin of `eps`.
 *
 * \def COMPARE_ENERGY
 * \ingroup test
 * Reduces the energy from the energy buffer and compares to the reference value
 * with a margin of `eps`.
 *
 * \def COMPARE_COUNT
 * \ingroup test
 * Reduces the number of interactions from the count buffer and compares to the
 * reference value.
 */
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


/**
 * \def COMPARE_VIR9
 * \ingroup test
 * Compares a virial tensor to the reference values with a margin of `eps`.
 *
 * \def COMPARE_VIR
 * \ingroup test
 * Reduces a virial tensor from a virial buffer and compares to the reference
 * values with a margin of `eps`.
 *
 * \def COMPARE_VIR2
 * \ingroup test
 * Reduces two virial tensors from two virial buffers and compares the sum to
 * the reference with a margin of `eps`.
 */
#define COMPARE_VIR9(vir1, ref_v, eps)                                         \
   {                                                                           \
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


/**
 * \def COMPARE_GRADIENT
 * \ingroup test
 * Copies out gradients from device to host and compares to the reference values
 * with a margin of `eps`.
 *
 * \def COMPARE_GRADIENT2
 * \ingroup test
 * Compares the flitered gradients[i][j] components.
 */
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


/**
 * \def COMPARE_BONDED_FORCE
 * \ingroup test
 * Compares the result of bonded (valence) force term.
 */
#define COMPARE_BONDED_FORCE(cpu_count, gpu_e, gpu_v, ref_e, eps_e, ref_count, \
                             ref_g, eps_g, ref_v, eps_v)                       \
   {                                                                           \
      auto do_ij_ = [](int, int) { return true; };                             \
      energy(calc::v3);                                                        \
      COMPARE_ENERGY(gpu_e, ref_e, eps_e);                                     \
      REQUIRE(cpu_count == ref_count);                                         \
                                                                               \
      energy(calc::v1);                                                        \
      COMPARE_ENERGY(gpu_e, ref_e, eps_e);                                     \
      COMPARE_GRADIENT2(ref_g, eps_g, do_ij_);                                 \
      COMPARE_VIR(gpu_v, ref_v, eps_v);                                        \
                                                                               \
      energy(calc::v4);                                                        \
      COMPARE_ENERGY(gpu_e, ref_e, eps_e);                                     \
      COMPARE_GRADIENT2(ref_g, eps_g, do_ij_);                                 \
                                                                               \
      energy(calc::v5);                                                        \
      COMPARE_GRADIENT2(ref_g, eps_g, do_ij_);                                 \
                                                                               \
      energy(calc::v6);                                                        \
      COMPARE_GRADIENT2(ref_g, eps_g, do_ij_);                                 \
      COMPARE_VIR(gpu_v, ref_v, eps_v);                                        \
   }
