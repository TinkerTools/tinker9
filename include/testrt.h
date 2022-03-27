#pragma once
#include "ff/atom.h"
#include "ff/energy.h"
#include "mod/md.h"
#include <string>
#include <vector>

namespace tinker {
/// \ingroup test
/// \brief Writes a file to disk in its constructor and removes this file in its destructor,
/// unless the file is set to be kept.
class TestFile
{
private:
   std::string name;
   bool good;

public:
   /// \brief Copies file from `src` to `dst` and append `extra` text to the file
   /// if `extra` is not empty. If `dst` is an empty string, the new file will
   /// be written to the current working directory with the same name as `src`.
   /// If `src` is an empty string, `dst` must be a valid name for the new file.
   TestFile(std::string src, std::string dst = "", std::string extra = "");
   /// \brief Removes the file on disk on exit.
   ~TestFile();
   /// \brief Prevents the file being deleted.
   void __keep();
};

/// \ingroup test
/// \brief Removes the file in its destructor as necessary.
class TestRemoveFileOnExit
{
private:
   std::string m_name;

public:
   TestRemoveFileOnExit(std::string fileToDelete);
   ~TestRemoveFileOnExit();
};

/// \ingroup test
/// \brief Reads reference values from a text file.
class TestReference
{
private:
   std::vector<double> gradient;
   double virial[3][3];
   double energy;
   int count;

public:
   TestReference(std::string pathToRefFile);
   int getCount() const;
   double getEnergy() const;
   const double (*getVirial() const)[3];
   const double (*getGradient() const)[3];
};

/// \ingroup test
/// \brief Returns tolerance eps depending on the predefined floating-point precision.
/// \param epsSingle  Larger `eps` for lower floating-point precision.
/// \param epsDouble  Smaller `eps` for higher floating-point precision.
double testGetEps(double epsSingle, double epsDouble);

/// \ingroup test
/// \brief Initializes the test.
void testBeginWithArgs(int argc, const char** argv);

/// \ingroup test
/// \brief Ends the test.
void testEnd();

/// \ingroup test
/// \brief Initializes MD in the test.
/// \param t    Temperature in Kelvin.
/// \param atm  Atmosphere in atm.
void testMdInit(double t = 0, double atm = 0);
}

/// \def COMPARE_INTS
/// \ingroup test
/// \brief Compare two integers.
///
/// \def COMPARE_REALS
/// \ingroup test
/// \brief Compare two floating-point numbers with a margin of `eps`.
///
/// \def COMPARE_ENERGY
/// \ingroup test
/// \brief Reduces the energy from the energy buffer and compares to the reference value
/// with a margin of `eps`.
///
/// \def COMPARE_COUNT
/// \ingroup test
/// \brief Reduces the number of interactions from the count buffer and compares to the
/// reference value.
#define COMPARE_INTS(i1, refi) REQUIRE(i1 == refi)
#define COMPARE_INTS_EPS(i1, refi, epsi)                                                           \
   {                                                                                               \
      int c1 = i1;                                                                                 \
      int r1 = refi;                                                                               \
      REQUIRE(r1 - epsi <= c1);                                                                    \
      REQUIRE(c1 <= r1 + epsi);                                                                    \
   }
#define COMPARE_REALS(v1, refv, eps) REQUIRE(v1 == Approx(refv).margin(eps))
#define COMPARE_ENERGY(gpuptr, ref_eng, eps)                                                       \
   {                                                                                               \
      double eng = energy_reduce(gpuptr);                                                          \
      REQUIRE(eng == Approx(ref_eng).margin(eps));                                                 \
   }
#define COMPARE_COUNT(gpuptr, ref_count)                                                           \
   {                                                                                               \
      int count = count_reduce(gpuptr);                                                            \
      REQUIRE(count == ref_count);                                                                 \
   }

/// \def COMPARE_VIR9
/// \ingroup test
/// \brief Compares a virial tensor to the reference values with a margin of `eps`.
///
/// \def COMPARE_VIR
/// \ingroup test
/// \brief Reduces a virial tensor from a virial buffer and compares to the reference
/// values with a margin of `eps`.
///
/// \def COMPARE_VIR2
/// \ingroup test
/// \brief Reduces two virial tensors from two virial buffers and compares the sum to
/// the reference with a margin of `eps`.
#define COMPARE_VIR9(vir1, ref_v, eps)                                                             \
   {                                                                                               \
      REQUIRE(vir1[0] == Approx(ref_v[0][0]).margin(eps));                                         \
      REQUIRE(vir1[1] == Approx(ref_v[0][1]).margin(eps));                                         \
      REQUIRE(vir1[2] == Approx(ref_v[0][2]).margin(eps));                                         \
      REQUIRE(vir1[3] == Approx(ref_v[1][0]).margin(eps));                                         \
      REQUIRE(vir1[4] == Approx(ref_v[1][1]).margin(eps));                                         \
      REQUIRE(vir1[5] == Approx(ref_v[1][2]).margin(eps));                                         \
      REQUIRE(vir1[6] == Approx(ref_v[2][0]).margin(eps));                                         \
      REQUIRE(vir1[7] == Approx(ref_v[2][1]).margin(eps));                                         \
      REQUIRE(vir1[8] == Approx(ref_v[2][2]).margin(eps));                                         \
   }
#define COMPARE_VIR(gpuptr, ref_v, eps)                                                            \
   {                                                                                               \
      virial_prec vir1[9];                                                                         \
      virial_reduce(vir1, gpuptr);                                                                 \
      REQUIRE(vir1[0] == Approx(ref_v[0][0]).margin(eps));                                         \
      REQUIRE(vir1[1] == Approx(ref_v[0][1]).margin(eps));                                         \
      REQUIRE(vir1[2] == Approx(ref_v[0][2]).margin(eps));                                         \
      REQUIRE(vir1[3] == Approx(ref_v[1][0]).margin(eps));                                         \
      REQUIRE(vir1[4] == Approx(ref_v[1][1]).margin(eps));                                         \
      REQUIRE(vir1[5] == Approx(ref_v[1][2]).margin(eps));                                         \
      REQUIRE(vir1[6] == Approx(ref_v[2][0]).margin(eps));                                         \
      REQUIRE(vir1[7] == Approx(ref_v[2][1]).margin(eps));                                         \
      REQUIRE(vir1[8] == Approx(ref_v[2][2]).margin(eps));                                         \
   }
#define COMPARE_VIR2(gpuptr, gpuptr2, ref_v, eps)                                                  \
   {                                                                                               \
      virial_prec vir1[9], vir2[9];                                                                \
      virial_reduce(vir1, gpuptr);                                                                 \
      virial_reduce(vir2, gpuptr2);                                                                \
      REQUIRE(vir1[0] + vir2[0] == Approx(ref_v[0][0]).margin(eps));                               \
      REQUIRE(vir1[1] + vir2[1] == Approx(ref_v[0][1]).margin(eps));                               \
      REQUIRE(vir1[2] + vir2[2] == Approx(ref_v[0][2]).margin(eps));                               \
      REQUIRE(vir1[3] + vir2[3] == Approx(ref_v[1][0]).margin(eps));                               \
      REQUIRE(vir1[4] + vir2[4] == Approx(ref_v[1][1]).margin(eps));                               \
      REQUIRE(vir1[5] + vir2[5] == Approx(ref_v[1][2]).margin(eps));                               \
      REQUIRE(vir1[6] + vir2[6] == Approx(ref_v[2][0]).margin(eps));                               \
      REQUIRE(vir1[7] + vir2[7] == Approx(ref_v[2][1]).margin(eps));                               \
      REQUIRE(vir1[8] + vir2[8] == Approx(ref_v[2][2]).margin(eps));                               \
   }

/// \def COMPARE_GRADIENT
/// \ingroup test
/// \brief Copies out gradients from device to host and compares to the reference values
/// with a margin of `eps`.
///
/// \def COMPARE_GRADIENT2
/// \ingroup test
/// \brief Compares the flitered gradients[i][j] components.
#define COMPARE_GRADIENT2(ref_grad, eps, check_ij)                                                 \
   {                                                                                               \
      std::vector<double> gradx(n), grady(n), gradz(n);                                            \
      copy_gradient(calc::grad, gradx.data(), grady.data(), gradz.data());                         \
      for (int i = 0; i < n; ++i) {                                                                \
         if (check_ij(i, 0))                                                                       \
            REQUIRE(gradx[i] == Approx(ref_grad[i][0]).margin(eps));                               \
         if (check_ij(i, 1))                                                                       \
            REQUIRE(grady[i] == Approx(ref_grad[i][1]).margin(eps));                               \
         if (check_ij(i, 2))                                                                       \
            REQUIRE(gradz[i] == Approx(ref_grad[i][2]).margin(eps));                               \
      }                                                                                            \
   }
#define COMPARE_GRADIENT(ref_grad, eps)                                                            \
   COMPARE_GRADIENT2(ref_grad, eps, [](int, int) { return true; })
