#include "mod/evalence.h"
#include "test/test.h"
#include "test/testrt.h"

using namespace tinker;

TEST_CASE("Urey-Ten-Water", "[ff][eurey][h2o10]")
{
   const char* k = "test_h2o10.key";
   const char* k0 = "ureyterm  only\n";
   const char* x1 = "test_h2o10.xyz";

   TestFile fke(TINKER9_DIRSTR "/src/test/file/water10/h2o10.key", k, k0);
   TestFile fx1(TINKER9_DIRSTR "/src/test/file/water10/h2o10.xyz", x1);
   TestFile fpr(TINKER9_DIRSTR "/src/test/file/commit_6fe8e913/water03.prm");

   TestReference r(TINKER9_DIRSTR "/src/test/ref/urey.txt");
   auto ref_e = r.get_energy();
   auto ref_v = r.get_virial();
   auto ref_count = r.get_count();
   auto ref_g = r.get_gradient();

   const char* argv[] = {"dummy", x1};
   int argc = 2;
   testBeginWithArgs(argc, argv);
   rc_flag = calc::xyz | calc::vmask;
   initialize();

   const double eps_e = 0.0001;
   const double eps_g = 0.0001;
   const double eps_v = 0.001;

   energy(calc::v3);
   COMPARE_REALS(esum, ref_e, eps_e);
   COMPARE_INTS(nurey, ref_count);

   energy(calc::v1);
   COMPARE_REALS(esum, ref_e, eps_e);
   COMPARE_GRADIENT(ref_g, eps_g);
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         COMPARE_REALS(vir[i * 3 + j], ref_v[i][j], eps_v);

   energy(calc::v4);
   COMPARE_REALS(esum, ref_e, eps_e);
   COMPARE_GRADIENT(ref_g, eps_g);

   energy(calc::v5);
   COMPARE_GRADIENT(ref_g, eps_g);

   energy(calc::v6);
   COMPARE_GRADIENT(ref_g, eps_g);
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         COMPARE_REALS(vir[i * 3 + j], ref_v[i][j], eps_v);

   finish();
   testEnd();
}
