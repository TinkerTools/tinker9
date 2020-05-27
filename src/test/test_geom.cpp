#include "files.h"
#include "test.h"
#include "test_rt.h"


using namespace tinker;


static const char* ortho_box = R"**(
a-axis 30.0
)**";


static const char* group_restrn = R"**(
group 1     5
group 2   7,8
group 3 -9,14

# force dist1 (dist2)
restrain-groups 1 2  2.1     1.3
restrain-groups 2 3  5.3     2.9 3.7
restrainterm        only
)**";


TEST_CASE("Geom-Local-Frame2-1", "[ff][egeom][local-frame2]")
{
   const char* k = "test_local_frame2.key";
   const char* x = "test_local_frame2.xyz";
   const char* argv[] = {"dummy", x};
   int argc = 2;
   int usage = calc::xyz | calc::mass | calc::vmask;


   std::string k0 = local_frame_key;
   k0 += ortho_box;


   TestFile fpr("amoeba09.prm", commit_6fe8e913::amoeba09_prm);
   TestFile fx1(x, local_frame_xyz2);


   SECTION("  - group restraint")
   {
      k0 += group_restrn;
      TestFile fke(k, k0);


      test_begin_with_args(argc, argv);
      rc_flag = usage;
      initialize();


      const double eps_e = 0.0001;
      const double ref_e = 11.0900;
      const int ref_count = 2;
      const double eps_g = 0.0001;
      const double eps_v = 0.001;
      const double ref_v[][3] = {{12.479, 12.738, 12.921},
                                 {12.738, 13.378, 16.537},
                                 {12.921, 16.537, 43.236}};
      const double ref_g[][3] = {
         {0.0000, 0.0000, 0.0000},    {0.0000, 0.0000, 0.0000},
         {0.0000, 0.0000, 0.0000},    {0.0000, 0.0000, 0.0000},
         {-5.0508, -4.7929, -1.9947}, {0.0000, 0.0000, 0.0000},
         {3.5223, 3.8316, 5.7537},    {3.5223, 3.8316, 5.7537},
         {-0.6803, -0.9793, -3.2457}, {-0.5833, -0.8398, -2.7832},
         {-0.5833, -0.8398, -2.7832}, {-0.0490, -0.0705, -0.2336},
         {-0.0490, -0.0705, -0.2336}, {-0.0490, -0.0705, -0.2336},
         {0.0000, 0.0000, 0.0000},    {0.0000, 0.0000, 0.0000},
         {0.0000, 0.0000, 0.0000},    {0.0000, 0.0000, 0.0000}};


      COMPARE_BONDED_FORCE(ngfix, eg, vir_eg, ref_e, eps_e, ref_count, ref_g,
                           eps_g, ref_v, eps_v);


      finish();
      test_end();
   }
}
