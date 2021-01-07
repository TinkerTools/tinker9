#include "test.h"
#include "test_rt.h"


using namespace tinker;


TEST_CASE("Improp-Trpcage", "[ff][eimprop][trpcage]")
{
   const char* k1 = "test_trpcage.key";
   const char* k2 = "impropterm only\n";
   const char* x1 = "test_trpcage.xyz";


   TestFile fke(TINKER9_DIRSTR "/src/test/file/trpcage/trp_charmm.key", k1, k2);
   TestFile fxy(TINKER9_DIRSTR "/src/test/file/trpcage/trp_charmm.xyz", x1);
   TestFile fpr(TINKER9_DIRSTR "/src/test/file/commit_11e84c69/charmm19.prm");


   TestReference r(TINKER9_DIRSTR "/src/test/ref/improp.txt");
   auto ref_e = r.get_energy();
   auto ref_v = r.get_virial();
   auto ref_count = r.get_count();
   auto ref_g = r.get_gradient();


   const double eps_e = test_get_eps(0.0001, 0.0001);
   const double eps_g = test_get_eps(0.0003, 0.0001);
   const double eps_v = test_get_eps(0.001, 0.001);


   const char* argv[] = {"dummy", x1};
   int argc = 2;
   test_begin_with_args(argc, argv);
   rc_flag = calc::xyz | calc::vmask;
   initialize();


   energy(calc::v3);
   COMPARE_REALS(esum, ref_e, eps_e);
   COMPARE_INTS(niprop, ref_count);


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
   test_end();
}
