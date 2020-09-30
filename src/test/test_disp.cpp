#include "test.h"
#include "test_rt.h"


using namespace tinker;


TEST_CASE("EDISP-1-NONDEWALD", "[ff][edisp][nondewald]")
{
   TestFile fx1(TINKER_GPU_DIRSTR "/src/test/file/c5h12acnh2.xyz");
   TestFile fk1(TINKER_GPU_DIRSTR "/src/test/file/disp/ndewald.key");
   TestFile fp1(TINKER_GPU_DIRSTR "/src/test/file/disp/hippo19.prm");
   const char* xn = "c5h12acnh2.xyz";
   const char* kn = "ndewald.key";
   const char* argv[] = {"dummy", xn, "-k", kn};
   int argc = 4;


   const double eps_e = 0.0001;
   const double eps_g = 0.0001;
   const double eps_v = 0.001;


   TestReference r(TINKER_GPU_DIRSTR "/src/test/disp.1.txt");
   auto ref_c = r.get_count();
   auto ref_e = r.get_energy();
   auto ref_v = r.get_virial();
   auto ref_g = r.get_gradient();


   rc_flag = calc::xyz | calc::vmask;
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


   energy(calc::v3);
   COMPARE_REALS(esum, ref_e, eps_e);
   COMPARE_INTS(count_reduce(ndisp), ref_c);


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


TEST_CASE("EDISP-2-DEWALD", "[ff][edisp][dewald]")
{
   TestFile fx1(TINKER_GPU_DIRSTR "/src/test/file/c5h12acnh2.xyz");
   TestFile fk1(TINKER_GPU_DIRSTR "/src/test/file/disp/dewald.key");
   TestFile fp1(TINKER_GPU_DIRSTR "/src/test/file/disp/hippo19.prm");
   const char* xn = "c5h12acnh2.xyz";
   const char* kn = "dewald.key";
   const char* argv[] = {"dummy", xn, "-k", kn};
   int argc = 4;


   const double eps_e = 0.001;
   const double eps_g = 0.0001;
   const double eps_v = 0.001;


   TestReference r(TINKER_GPU_DIRSTR "/src/test/disp.2.txt");
   auto ref_c = r.get_count();
   auto ref_e = r.get_energy();
   auto ref_v = r.get_virial();
   auto ref_g = r.get_gradient();


   rc_flag = calc::xyz | calc::vmask;
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


   energy(calc::v3);
   COMPARE_REALS(esum, ref_e, eps_e);
   COMPARE_INTS(count_reduce(ndisp), ref_c);


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
