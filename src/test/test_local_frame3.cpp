#include "files.h"
#include "test.h"
#include "test_rt.h"


using namespace TINKER_NAMESPACE;


namespace {
int usage = calc::xyz | calc::vmask;
const char* key0 = R"**(
parameters  amoeba09
cutoff           7.0
neighbor-list
list-buffer      0.1
a-axis          20.0


bondterm        none
angleterm       none
strbndterm      none
ureyterm        none
torsionterm     none
vdwterm         none
)**";
const char* kname = "test_local_frame3.key";
const char* xname = "test_local_frame3.xyz";
const char* argv[] = {"dummy", xname};
int argc = 2;
}


TEST_CASE("Local-Frame3-1",
          "[ff][empole][epolar][emplar][nonewald][local-frame3]")
{
   std::string key = key0;

   TestFile fxy(xname, local_frame_xyz2);
   TestFile fke(kname, key);
   TestFile fpr("amoeba09.prm", commit_6fe8e913::amoeba09_prm);

   test_begin_with_args(argc, argv);
   rc_flag = usage;
   initialize();

   SECTION("  - emplar -- non-ewald pbc")
   {
      const int ref_count = 120;
      // empole -107.5885 count 120
      // epolar  -18.4063 count 120
      const double ref_eng = -125.9948;
      // empole virial
      //   83.083   6.570 -33.345
      //    6.570   2.497 -13.213
      //  -33.345 -13.213  38.666
      // epolar virial
      //   31.264  -0.788   7.542
      //   -0.788   4.922  -3.701
      //    7.542  -3.701  18.821
      const double ref_v[][3] = {{114.347, 5.782, -25.803},
                                 {5.782, 7.419, -16.914},
                                 {-25.803, -16.914, 57.488}};
      const double ref_g[][3] = {
         {21.9813, 3.0377, -6.1168},  {3.2451, -11.5492, 18.3287},
         {9.1717, 3.2256, -2.2598},   {-3.8440, -1.5729, 0.9964},
         {-4.1469, -2.1011, 0.0910},  {16.3362, -14.6835, 3.6483},
         {-5.9219, 5.8674, 3.6770},   {-25.8555, 11.8977, -4.4568},
         {1.1900, -0.5955, 9.1827},   {-0.5987, -0.2092, -3.9256},
         {3.4362, -0.4291, 2.6352},   {-1.6781, 0.2599, -1.3334},
         {-2.3079, 0.2314, -1.4129},  {-2.0652, 0.8646, -1.6531},
         {-28.0256, -5.4435, 3.8673}, {2.0955, 2.9518, -1.4402},
         {8.6696, -0.2665, -2.3909},  {8.3182, 8.5144, -17.4371}};
      const double eps_e = 0.0001;
      const double eps_g = 0.0001;
      const double eps_v = 0.001;

      energy(calc::v0);
      COMPARE_REALS(esum, ref_eng, eps_e);

      energy(calc::v1);
      COMPARE_REALS(esum, ref_eng, eps_e);
      COMPARE_GRADIENT_(ref_g, eps_g);
      for (int i = 0; i < 3; ++i)
         for (int j = 0; j < 3; ++j)
            COMPARE_REALS(vir[i * 3 + j], ref_v[i][j], eps_v);

      energy(calc::v3);
      COMPARE_REALS(esum, ref_eng, eps_e);
      COMPARE_INTS(get_count(nem), ref_count);
      COMPARE_INTS(get_count(nep), ref_count);

      energy(calc::v4);
      COMPARE_REALS(esum, ref_eng, eps_e);
      COMPARE_GRADIENT_(ref_g, eps_g);

      energy(calc::v5);
      COMPARE_GRADIENT_(ref_g, eps_g);

      energy(calc::v6);
      COMPARE_GRADIENT_(ref_g, eps_g);
      for (int i = 0; i < 3; ++i)
         for (int j = 0; j < 3; ++j)
            COMPARE_REALS(vir[i * 3 + j], ref_v[i][j], eps_v);
   }

   finish();
   test_end();
}


TEST_CASE("Local-Frame3-2", "[ff][empole][epolar][emplar][ewald][local-frame3]")
{
   std::string key = key0;
   key += "ewald\n";

   TestFile fxy(xname, local_frame_xyz2);
   TestFile fke(kname, key);
   TestFile fpr("amoeba09.prm", commit_6fe8e913::amoeba09_prm);

   test_begin_with_args(argc, argv);
   rc_flag = usage;
   initialize();

   SECTION("  - emplar -- ewald pbc") {}

   finish();
   test_end();
}
