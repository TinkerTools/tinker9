#include "test.h"
#include "testrt.h"
#include "tinker9.h"
#include "tool/ioprint.h"

#include <tinker/detail/units.hh>
#include <tinker/routines.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

using namespace tinker;

namespace {
/// Ends the Fortran runtime however the case ends.
struct FortranRuntimeScope
{
   FortranRuntimeScope(int argc, const char** argv)
   {
      tinkerFortranRuntimeBegin(argc, const_cast<char**>(argv));
      tinker_f_command();
   }
   ~FortranRuntimeScope() { tinkerFortranRuntimeEnd(); }
};
}

// CUDA builds only. The bar driver does not depend on the platform, and the
// OpenACC build would either rerun the CUDA kernels, its default for a key with
// a neighbor list and a box, or, under GPU_PACKAGE=OPENACC, spend minutes per
// case on the 12501 atom box.
#if TINKER_GPULANG_CUDA

namespace {
/// One trajectory block of a BAR file: its header, then one row per frame.
struct BarBlock
{
   int nfrm = 0;
   double temp = 0;
   std::vector<double> u0, u1, vol; ///< Energies in states 0 and 1, and the box volume (0 if none).
};

/// Reads the two trajectory blocks "bar" mode 1 writes.
std::vector<BarBlock> readBar(const std::string& path)
{
   std::ifstream fin(path);
   REQUIRE(fin);
   std::vector<BarBlock> blocks;
   std::string line;
   while (std::getline(fin, line)) {
      BarBlock b;
      std::istringstream hdr(line);
      if (not(hdr >> b.nfrm >> b.temp))
         continue;
      for (int i = 0; i < b.nfrm; ++i) {
         REQUIRE(std::getline(fin, line));
         std::istringstream row(line);
         int idx;
         double u0, u1, vol = 0;
         REQUIRE((row >> idx >> u0 >> u1));
         row >> vol;
         REQUIRE(idx == i + 1);
         b.u0.push_back(u0);
         b.u1.push_back(u1);
         b.vol.push_back(vol);
      }
      blocks.push_back(b);
   }
   return blocks;
}

// Runs "bar 1 <A>.arc 298 <B>.arc 298 N", which re-evaluates both trajectories
// in both states and writes <A>.bar, and compares that file with the Tinker
// reference (test_bar.f). Each AMOEBA water box trajectory decouples one water
// molecule between two neighboring lambda windows.
void runBar(const std::string& basea, const std::string& baseb)
{
   TestLmdaFlagReset lmdaReset;

   const std::string files = TINKER9_DIRSTR "/test/file/bar/";
   const std::string arca = basea + ".arc", arcb = baseb + ".arc";
   TestFile fprm(files + "hostsG4.prm");
   TestFile fkeya(files + basea + ".key");
   TestFile fkeyb(files + baseb + ".key");
   TestFile farca(files + arca);
   TestFile farcb(files + arcb);

   // A leftover output would make "bar" write <A>.bar_2 instead.
   const std::string barfile = basea + ".bar";
   std::remove(barfile.c_str());
   TestRemoveFileOnExit rmbar(barfile);

   const char* argv[] = {"dummy", "1", arca.c_str(), "298", arcb.c_str(), "298", "N"};
   const int argc = 7;
   {
      FortranRuntimeScope scope(argc, argv);
      xBar(argc, const_cast<char**>(argv));
   }

   auto got = readBar(barfile);
   auto ref = readBar(std::string(TINKER9_DIRSTR "/test/ref/bar/") + barfile);
   REQUIRE(got.size() == 2);
   REQUIRE(ref.size() == 2);

   // Total energies of 12501 atoms near -37000 kcal/mol. Catch2's Approx keeps
   // a relative slack of about 0.4 kcal/mol at that size, more than the state 1
   // minus state 0 differences BAR works from, so the comparisons are absolute.
   // Mixed precision lands about 0.02 from the reference in total energy and
   // 2e-4 in the difference, and the volume, computed in single precision, up
   // to 0.01 away; double precision agrees to the printed 1e-4. The difference
   // is the sensitive check: moving VDW-LAMBDA by 0.001 shifts it by 3e-3 or
   // more on every frame.
   const double eps_e = testGetEps(1.0e-1, 1.0e-3);
   const double eps_d = testGetEps(2.0e-3, 5.0e-4);
   const double eps_vol = testGetEps(5.0e-2, 1.0e-3);
   for (int k = 0; k < 2; ++k) {
      const auto &g = got[k], &r = ref[k];
      REQUIRE(g.nfrm == r.nfrm);
      REQUIRE(g.temp == r.temp);
      for (int i = 0; i < r.nfrm; ++i) {
         INFO("block " << k + 1 << " frame " << i + 1);
         REQUIRE(std::fabs(g.u0[i] - r.u0[i]) <= eps_e);
         REQUIRE(std::fabs(g.u1[i] - r.u1[i]) <= eps_e);
         REQUIRE(std::fabs((g.u1[i] - g.u0[i]) - (r.u1[i] - r.u0[i])) <= eps_d);
         REQUIRE(std::fabs(g.vol[i] - r.vol[i]) <= eps_vol);
      }
   }
}
}

TEST_CASE("BAR-water-00-0000-0950", "[ff][bar]") { runBar("water-00-0000-0950", "water-00-0000-0900"); }
TEST_CASE("BAR-water-00-0950-1000", "[ff][bar]") { runBar("water-00-0950-1000", "water-00-0900-1000"); }

#endif

// "bar" mode 2 reads a BAR file and evaluates no energies, so the cases below
// run on every platform. The bootstrap draws from a clock seed, so each case
// checks a property that holds whatever the seed.
namespace {
/// Sends stdout to a file until destroyed.
class StdoutToFile
{
   int m_saved;

public:
   StdoutToFile(const std::string& file)
   {
      std::fflush(stdout);
      m_saved = dup(fileno(stdout));
      REQUIRE(m_saved != -1);
      REQUIRE(std::freopen(file.c_str(), "w", stdout));
   }
   ~StdoutToFile()
   {
      std::fflush(stdout);
      dup2(m_saved, fileno(stdout));
      close(m_saved);
   }
};

/// Runs "bar 2 <file> 1 <na> 1 1 <nb> 1" over every frame and returns the
/// output.
std::string runBarCalc(const std::string& barfile, int nfrma, int nfrmb)
{
   const std::string logfile = "bar-mode2.log";
   TestRemoveFileOnExit rmlog(logfile);
   const std::string na = std::to_string(nfrma), nb = std::to_string(nfrmb);
   const char* argv[] = {"dummy", "2", barfile.c_str(), "1", na.c_str(), "1", "1", nb.c_str(), "1"};
   const int argc = 9;
   {
      StdoutToFile redirect(logfile);
      FortranRuntimeScope scope(argc, argv);
      xBar(argc, const_cast<char**>(argv));
   }

   std::ifstream fin(logfile);
   std::stringstream ss;
   ss << fin.rdbuf();
   return ss.str();
}

/// Reads "<value> +/- <sd>" after the first occurrence of label.
bool readEstimate(const std::string& out, const std::string& label, double& val, double& sd)
{
   auto pos = out.find(label);
   if (pos == std::string::npos)
      return false;
   return std::sscanf(out.c_str() + pos + label.size(), "%lf +/-%lf", &val, &sd) == 2;
}

/// One trajectory block at 298 K, with more digits than mode 1 writes. Rows
/// have no volume if vol is empty.
std::string barBlock(const std::vector<double>& u0, const std::vector<double>& u1, const std::vector<double>& vol)
{
   std::string s = format("%8d%10.2lf  Synthetic\n", (int)u0.size(), 298.0);
   for (size_t i = 0; i < u0.size(); ++i) {
      s += format("%8d  %18.10lf%18.10lf", (int)i + 1, u0[i], u1[i]);
      if (vol.size())
         s += format("%18.10lf", vol[i]);
      s += "\n";
   }
   return s;
}

const double rt298 = units::gasconst * 298.0;
const char* const barIter = " Free Energy via BAR Iteration";
const char* const barBoot = " Free Energy via BAR Bootstrap";
}

// The two BAR files "bar" mode 1 writes in the cases above. The BAR iteration
// is deterministic; the values are Tinker's bar.x mode 2 output.
TEST_CASE("BAR-mode2-ref", "[bar]")
{
   const struct
   {
      const char* base;
      double val, sd;
   } refs[] = {{"water-00-0000-0950", -0.2070, 0.4380}, {"water-00-0950-1000", 1.3773, 0.4279}};
   for (const auto& r : refs) {
      INFO(r.base);
      const std::string barfile = std::string(r.base) + ".bar";
      TestFile fbar(TINKER9_DIRSTR "/test/ref/bar/" + barfile);
      auto out = runBarCalc(barfile, 2, 2);
      double val, sd;
      REQUIRE(readEstimate(out, barIter, val, sd));
      REQUIRE(std::fabs(val - r.val) <= 1.0e-4);
      REQUIRE(std::fabs(sd - r.sd) <= 1.0e-4);
   }
}

// Every frame has dU = 0.5 in both directions, so BAR gives 0.5 exactly. Rows
// with no volume and rows with an explicit zero volume must both mean NVT.
TEST_CASE("BAR-mode2-zero-volume", "[bar]")
{
   const int n = 10;
   std::vector<double> ua0, ua1, ub0, ub1;
   for (int i = 0; i < n; ++i) {
      ua0.push_back(-100.0 + 0.37 * i);
      ua1.push_back(ua0[i] + 0.5);
      ub0.push_back(-90.0 + 0.21 * i);
      ub1.push_back(ub0[i] + 0.5);
   }
   const std::vector<double> none, zero(n, 0.0);
   for (const auto* vol : {&none, &zero}) {
      INFO((vol->size() ? "explicit zero volume" : "no volume column"));
      const std::string barfile = "bar-mode2-nvt.bar";
      TestFile fbar("", barfile, barBlock(ua0, ua1, *vol) + barBlock(ub0, ub1, *vol));
      auto out = runBarCalc(barfile, n, n);
      double val, sd;
      REQUIRE(readEstimate(out, barIter, val, sd));
      REQUIRE(std::fabs(val - 0.5) <= 1.0e-4);
      REQUIRE(out.find("PdV") == std::string::npos);
   }
}

// The volume correction is -RT ln(V/<V>). Making dU = c + RT ln(V) frame by
// frame leaves every frame's corrected dU equal to c + RT ln(<V>), so every
// bootstrap resample gives the same BAR value and the bootstrap deviation is
// about zero. That holds only if each frame's energies are corrected with its
// own volume; correcting with another frame's volume leaves a spread of up to
// RT ln 3 and a deviation of a few tenths.
TEST_CASE("BAR-mode2-bootstrap-volume-index", "[bar]")
{
   const int n = 10;
   std::vector<double> ua0, ua1, vola, ub0, ub1, volb;
   for (int i = 0; i < n; ++i) {
      vola.push_back(1000.0 + 2000.0 * i / (n - 1));
      ua0.push_back(-100.0 - 0.3 * i);
      ua1.push_back(ua0[i] + 0.5 + rt298 * std::log(vola[i] / 2000.0));
      volb.push_back(3000.0 - 1500.0 * i / (n - 1));
      ub1.push_back(-95.0 + 0.2 * i);
      ub0.push_back(ub1[i] - 0.5 + rt298 * std::log(volb[i] / 2000.0));
   }
   const std::string barfile = "bar-mode2-npt.bar";
   TestFile fbar("", barfile, barBlock(ua0, ua1, vola) + barBlock(ub0, ub1, volb));
   auto out = runBarCalc(barfile, n, n);
   double val, sd, bval, bsd;
   REQUIRE(readEstimate(out, barIter, val, sd));
   REQUIRE(readEstimate(out, barBoot, bval, bsd));
   REQUIRE(std::fabs(bval - val) <= 1.0e-3);
   // Rounding can make the variance a tiny negative number and print nan.
   REQUIRE_FALSE(bsd > 1.0e-3);
}

// Nine of the ten frames of A have no overlap with B, so a resample of A that
// misses the remaining frame drives the BAR iteration to inf and then NaN.
// About a third of the trials do, and they must be skipped, not loop forever.
TEST_CASE("BAR-mode2-bootstrap-nonconverged", "[bar]")
{
   const int n = 10;
   std::vector<double> ua0, ua1, ub0, ub1;
   for (int i = 0; i < n; ++i) {
      ua0.push_back(-100.0 - 0.1 * i);
      ua1.push_back(ua0[i] + (i == 0 ? 0.0 : 1000.0));
      ub1.push_back(-100.0 + 0.1 * i);
      ub0.push_back(ub1[i] + 0.3 * (i % 3) - 0.3);
   }
   const std::string barfile = "bar-mode2-fail.bar";
   TestFile fbar("", barfile, barBlock(ua0, ua1, {}) + barBlock(ub0, ub1, {}));
   auto out = runBarCalc(barfile, n, n);
   double val, sd, bval, bsd;
   REQUIRE(readEstimate(out, barIter, val, sd));
   const std::string warn = " BAR Bootstrap:  ";
   auto pos = out.find(warn);
   REQUIRE(pos != std::string::npos);
   int nfail = 0, ntrial = 0;
   REQUIRE(std::sscanf(out.c_str() + pos + warn.size(), "%d of %d", &nfail, &ntrial) == 2);
   REQUIRE(nfail > 0);
   REQUIRE(nfail < ntrial);
   REQUIRE(readEstimate(out, barBoot, bval, bsd));
   REQUIRE(std::isfinite(bval));
   REQUIRE(std::isfinite(bsd));
}
