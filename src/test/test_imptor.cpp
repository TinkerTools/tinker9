#include "test/test.h"
#include "test/testrt.h"

using namespace tinker;

TEST_CASE("ALA-1", "[ff][eimptor][ala]")
{
   const char* k = "test_ala.key";
   const char* k0 = "imptorterm   only\n"
                    "imptorunit    100\n";
   const char* x = "test_ala.xyz";

   TestFile fpr(TINKER9_DIRSTR "/src/test/file/commit_350df099/amber99sb.prm");
   TestFile fx1(TINKER9_DIRSTR "/src/test/file/ala/ala.xyz", x);
   TestFile fk1(TINKER9_DIRSTR "/src/test/file/ala/ala.key", k, k0);

   TestReference r(TINKER9_DIRSTR "/src/test/ref/imptor.txt");
   auto ref_e = r.get_energy();
   auto ref_v = r.get_virial();
   auto ref_count = r.get_count();
   auto ref_g = r.get_gradient();

   const double eps_e = testGetEps(0.0002, 0.0001);
   const double eps_g = 0.0001;
   const double eps_v = 0.001;

   const char* argv[] = {"dummy", x};
   int argc = 2;
   int usage = calc::xyz | calc::vmask;
   testBeginWithArgs(argc, argv);
   rc_flag = usage;
   initialize();

   energy(calc::v3);
   COMPARE_REALS(esum, ref_e, eps_e);
   COMPARE_INTS(nitors, ref_count);

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
