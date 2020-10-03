#include "files.h"
#include "test.h"
#include "test_rt.h"


using namespace tinker;


TEST_CASE("ALA-1", "[ff][eimptor][ala]")
{
   const char* k = "test_ala.key";
   const char* x = "test_ala.xyz";
   const char* argv[] = {"dummy", x};
   int argc = 2;
   int usage = calc::xyz | calc::vmask;

   std::string k0 = ala_amber99sb_key;
   k0 += "imptorterm   only\n"
         "imptorunit    100\n";


   TestFile fpr("amber99sb.prm", commit_350df099::amber99sb_prm);
   TestFile fx1(x, ala_amber99sb_xyz);
   TestFile fk1(k, k0);


   test_begin_with_args(argc, argv);
   rc_flag = usage;
   initialize();


   const double eps_e = 0.0001;
   const double ref_e = 1.64524;
   const int ref_count = 4;
   const double eps_g = 0.0001;
   const double eps_v = 0.001;
   const double ref_v[][3] = {{-2.1348, -2.0662, -0.6327},
                              {-2.0662, 2.1204, 0.0670},
                              {-0.6327, 0.0670, 0.0144}};
   const double ref_g[][3] = {
      {1.1233, -4.6732, 0.0347},       {2.0893, -10.0294, -0.5275},
      {1.4803, -6.1775, 0.0372},       {0.0000, 0.0000, 0.0000},
      {0.0000, 0.0000, 0.0000},        {0.0000, 0.0000, 0.0000},
      {-20.3880, 75.5182, 0.9391},     {-24.3578, 12.1263, 14.2926},
      {104.7433, -124.1869, -51.0335}, {-39.4244, 47.4284, 19.4670},
      {9.6264, -31.3594, -0.2040},     {0.0000, 0.0000, 0.0000},
      {0.0000, 0.0000, 0.0000},        {0.0000, 0.0000, 0.0000},
      {0.0000, 0.0000, 0.0000},        {0.0000, 0.0000, 0.0000},
      {-33.8757, 40.1384, 16.4978},    {-0.4310, 0.5154, 0.2110},
      {-0.5858, 0.6998, 0.2856},       {0.0000, 0.0000, 0.0000},
      {0.0000, 0.0000, 0.0000},        {0.0000, 0.0000, 0.0000}};


   COMPARE_BONDED_FORCE(nitors, eit, vir_eit, ref_e, eps_e, ref_count, ref_g,
                        eps_g, ref_v, eps_v);


   finish();
   test_end();
}
