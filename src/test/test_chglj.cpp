#include "test.h"
#include "test_rt.h"


using namespace tinker;


TEST_CASE("Chglj-Trpcage", "[ff][echarge][evdw][echglj][lj][trpcage]")
{
   rc_flag = calc::xyz | calc::energy | calc::grad | calc::virial;


   const char* kn = "test_chglj.key";
   std::string k0 = R"**(
bondterm     none
angleterm    none
impropterm   none
torsionterm  none
)**";
   const char* xn = "test_chglj.xyz";


   const double eps_e = 0.0001;
   const double eps_g = 0.0005;
   const double eps_v = 0.001;
   const char* argv[] = {"dummy", xn};
   int argc = 2;


   SECTION("  - echglj -- no pbc, no cutoff")
   {
      TestFil2 fx1(TINKER9_DIRSTR "/src/test/file/trpcage/trp_charmm.xyz", xn);
      TestFil2 fk1(TINKER9_DIRSTR "/src/test/file/trpcage/trp_charmm.key", kn,
                   k0);
      TestFile fp1(TINKER9_DIRSTR
                   "/src/test/file/commit_11e84c69/charmm19.prm");


      TestReference r(TINKER9_DIRSTR "/src/test/chglj.1.txt");
      auto ref_e = r.get_energy();
      auto ref_v = r.get_virial();
      auto ref_g = r.get_gradient();


      test_begin_with_args(argc, argv);
      initialize();


      energy(calc::v0);
      COMPARE_REALS(esum, ref_e, eps_e);


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


   SECTION("  - echglj -- pbc, cutoff")
   {
      std::string k1 = k0 + R"**(
ewald
cutoff 9.0
neighbor-list
list-buffer 0.5
a-axis  30
b-axis  25
c-axis  20
vdw-correction
)**";


      TestFil2 fx1(TINKER9_DIRSTR "/src/test/file/trpcage/trp_charmm.xyz", xn);
      TestFil2 fk1(TINKER9_DIRSTR "/src/test/file/trpcage/trp_charmm.key", kn,
                   k1);
      TestFile fp1(TINKER9_DIRSTR
                   "/src/test/file/commit_11e84c69/charmm19.prm");


      TestReference r(TINKER9_DIRSTR "/src/test/chglj.2.txt");
      auto ref_e = r.get_energy();
      auto ref_v = r.get_virial();
      auto ref_g = r.get_gradient();


      test_begin_with_args(argc, argv);
      initialize();


      energy(calc::v0);
      COMPARE_REALS(esum, ref_e, eps_e);


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
}
