#include "test/test.h"
#include "test/testrt.h"

using namespace tinker;

TEST_CASE("ECHGTRN-1", "[ff][hippo][echgtrn][dmso]")
{
   TestFile fx1(TINKER9_DIRSTR "/src/test/file/dmso/dmso.xyz");
   TestFile fk1(TINKER9_DIRSTR "/src/test/file/hippo/chgtrn/chgtrn.key");
   TestFile fp1(TINKER9_DIRSTR "/src/test/file/hippo/hippo19.prm");
   const char* xn = "dmso.xyz";
   const char* kn = "chgtrn.key";
   const char* argv[] = {"dummy", xn, "-k", kn};
   int argc = 4;

   const double eps_e = testGetEps(0.0001, 0.0001);
   const double eps_g = testGetEps(0.0001, 0.0001);
   const double eps_v = testGetEps(0.0005, 0.0004);

   TestReference r(TINKER9_DIRSTR "/src/test/ref/chgtrn.1.txt");
   auto ref_c = r.get_count();
   auto ref_e = r.get_energy();
   auto ref_v = r.get_virial();
   auto ref_g = r.get_gradient();

   rc_flag = calc::xyz | calc::vmask;
   testBeginWithArgs(argc, argv);
   initialize();

   energy(calc::v0);
   COMPARE_REALS(esum, ref_e, eps_e);

   energy(calc::v1);
   COMPARE_REALS(esum, ref_e, eps_e);
   COMPARE_GRADIENT(ref_g, eps_g);
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         COMPARE_REALS(vir[i * 3 + j], ref_v[i][j], eps_v);

   energy(calc::v3);
   COMPARE_REALS(esum, ref_e, eps_e);
   COMPARE_INTS(count_reduce(nct), ref_c);

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
