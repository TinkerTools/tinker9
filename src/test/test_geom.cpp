#include "test.h"
#include "test_rt.h"


using namespace tinker;


TEST_CASE("Geom-Local-Frame2-1", "[ff][egeom][local-frame2]")
{
   const char* k = "test_local_frame2.key";
   const char* x = "test_local_frame2.xyz";
   const char* argv[] = {"dummy", x};
   int argc = 2;
   int usage = calc::xyz | calc::mass | calc::vmask;


   TestFile fpr(TINKER9_DIRSTR "/src/test/file/commit_6fe8e913/amoeba09.prm");
   TestFile fx1(TINKER9_DIRSTR "/src/test/file/local_frame/local_frame2.xyz",
                x);


   SECTION("  - group restraint")
   {
      const char* k0 =
         R"**(
a-axis   30.0

group 1     5
group 2   7,8
group 3 -9,14

# force dist1 (dist2)
restrain-groups 1 2  2.1     1.3
restrain-groups 2 3  5.3     2.9 3.7
restrainterm        only
)**";
      TestFile fke(TINKER9_DIRSTR "/src/test/file/local_frame/local_frame.key",
                   k, k0);


      TestReference r(TINKER9_DIRSTR "/src/test/ref/geom.1.txt");
      auto ref_e = r.get_energy();
      auto ref_v = r.get_virial();
      auto ref_count = r.get_count();
      auto ref_g = r.get_gradient();


      const double eps_e = 0.0001;
      const double eps_g = 0.0001;
      const double eps_v = 0.001;


      test_begin_with_args(argc, argv);
      rc_flag = usage;
      initialize();


      energy(calc::v3);
      COMPARE_REALS(esum, ref_e, eps_e);
      COMPARE_INTS(ngfix, ref_count);


      energy(calc::v1);
      COMPARE_REALS(esum, ref_e, eps_e);
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
