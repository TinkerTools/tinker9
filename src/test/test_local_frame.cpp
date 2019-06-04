#include "files.h"
#include "test/ff.h"
#include "test/rt.h"
#include "test/test.h"

m_tinker_using_namespace;
using namespace test;

TEST_CASE("Local-Frame-1", "[ff][empole][coulomb][local-frame]") {
  file fpr("amoeba09.prm", amoeba09_prm);

  const char* k = "test_local_frame.key";
  std::string key0 = local_frame_key;

  const char* x1 = "test_local_frame.xyz";
  file fx1(x1, local_frame_xyz);

  int usage = 0;
  usage |= gpu::use_xyz;
  usage |= gpu::use_energy;
  usage |= gpu::use_grad;
  usage |= gpu::use_virial;

  SECTION("empole -- gas phase, no cutoff") {
    std::string key1 = key0;
    key1 += "multipoleterm    only\n";
    file fke(k, key1);

    const char* argv[] = {"dummy", x1};
    int argc = 2;

    const double eps_e = 1.0e-5;
    const double ref_eng = -95.7432;
    const int ref_count = 201;
    const double eps_g = 0.0005;
    const double ref_grad[][3] = {
        {-20.9915, 29.4368, 0.5977},   {12.7503, -9.9272, 3.3604},
        {15.0932, -5.2916, 0.1809},    {-6.0154, 0.6960, -0.2064},
        {-5.0594, 4.3698, -1.3173},    {36.4924, -26.8619, 11.0059},
        {-11.2357, 10.4015, -1.8246},  {-26.9426, 13.9832, -5.6296},
        {12.3270, -4.4632, 9.1612},    {-6.1247, 3.8477, -1.6673},
        {3.8365, -3.0201, -1.0674},    {-1.5792, 0.9991, -0.5125},
        {-2.3155, 1.5935, 0.3036},     {-2.0629, 1.7467, 0.9395},
        {0.1139, -5.1353, 9.9238},     {-1.9284, 1.1291, -2.8500},
        {1.1353, -0.2388, -5.6115},    {-0.2273, 6.5552, -7.8322},
        {23.9256, -41.7121, -54.8522}, {-18.4759, 10.3134, 8.4037},
        {-6.0110, 9.7698, 14.0093},    {3.2955, 1.8082, 25.4850},
    };
    const double eps_v = 0.005;
    const double ref_v[][3] = {{72.636, -34.992, 9.243},
                               {-34.992, 53.410, 6.517},
                               {9.243, 6.517, 4.708}};

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    gpu::zero_egv();
    gpu::elec_init(gpu::v0);
    tinker_gpu_empole0();
    gpu::torque(gpu::v0);
    COMPARE_ENERGY_(gpu::em, ref_eng, eps_e);

    gpu::zero_egv();
    gpu::elec_init(gpu::v1);
    tinker_gpu_empole1();
    gpu::torque(gpu::v1);
    COMPARE_ENERGY_(gpu::em, ref_eng, eps_e);
    COMPARE_GRADIENT_(ref_grad, eps_g);
    COMPARE_VIR2_(gpu::vir_em, gpu::vir_trq, ref_v, eps_v);

    gpu::zero_egv();
    tinker_gpu_empole3();
    COMPARE_ENERGY_(gpu::em, ref_eng, eps_e);
    COMPARE_COUNT_(gpu::nem, ref_count);

    tinker_gpu_data_destroy();
    test_end();
  }
}

TEST_CASE("Local-Frame-2", "[ff][empole][ewald][local-frame]") {
  file fpr("amoeba09.prm", amoeba09_prm);

  const char* k = "test_local_frame.key";
  std::string key0 = local_frame_key;

  const char* x1 = "test_local_frame.xyz";
  file fx1(x1, local_frame_xyz);

  int usage = 0;
  usage |= gpu::use_xyz;
  usage |= gpu::use_energy;
  usage |= gpu::use_grad;
  usage |= gpu::use_virial;

  SECTION("empole -- pme") {
    std::string key1 = key0;
    key1 += "multipoleterm    only\n";
    key1 += "ewald\n";
    key1 += "ewald-cutoff    7.0\n";
    key1 += "neighbor-list\n";
    key1 += "list-buffer    0.1\n";
    key1 += "a-axis    20.0\n";
    file fke(k, key1);

    const char* argv[] = {"dummy", x1};
    int argc = 2;

    const double eps_e = 0.0001;
    // emreal 132.76851908903473
    // +recip 279.93038602122442 // recip 147.1618669321897
    // +self -99.229836789909172 // self -379.1602228111336
    const double ref_ereal = 132.76851908903473;
    const double ref_erecip = 147.1618669321897;
    const double ref_eself = -379.1602228111336;
    const double ref_eng = ref_eself + ref_erecip + ref_ereal;
    const int ref_count = 222;
    const double eps_g = 0.0005;
    const double ref_grad[][3] = {
        {-20.0239, 29.1238, 0.7844},   {11.8249, -9.4636, 3.1958},
        {14.5683, -4.9576, 0.0931},    {-5.7229, 0.5038, -0.1740},
        {-4.7824, 4.2038, -1.2616},    {35.8953, -26.6166, 10.8531},
        {-10.9482, 10.2634, -1.7619},  {-26.6406, 13.8783, -5.5468},
        {11.7934, -4.2057, 8.9439},    {-5.7126, 3.6221, -1.4949},
        {3.5447, -2.8759, -1.0921},    {-1.4544, 0.9341, -0.4978},
        {-2.1774, 1.4935, 0.3237},     {-1.9210, 1.7174, 0.9565},
        {-0.6072, -4.9773, 9.5574},    {-1.8037, 1.1493, -2.7232},
        {1.4169, -0.3143, -5.4280},    {0.0561, 6.4388, -7.8332},
        {22.8962, -41.3767, -55.0242}, {-17.8085, 9.9950, 8.4640},
        {-5.9679, 9.6342, 14.0679},    {3.6422, 2.0227, 25.5414}};
    const double eps_v = 0.001;
    const double ref_v[][3] = {{76.909, -31.146, 7.610},
                               {-31.146, 53.934, 7.245},
                               {7.610, 7.245, 3.181}};

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    gpu::zero_egv();
    gpu::elec_init(gpu::v1);
    tinker_gpu_empole1();
    gpu::torque(gpu::v1);
    COMPARE_ENERGY_(gpu::em, ref_eng, eps_e);
    COMPARE_GRADIENT_(ref_grad, eps_g);
    COMPARE_VIR2_(gpu::vir_em, gpu::vir_trq, ref_v, eps_v);

    tinker_gpu_data_destroy();
    test_end();
  }
}

TEST_CASE("Local-Frame-3", "[ff][epolar][ewald][local-frame]") {
  file fpr("amoeba09.prm", amoeba09_prm);

  const char* k = "test_local_frame.key";
  std::string key0 = local_frame_key;

  const char* x1 = "test_local_frame.xyz";
  file fx1(x1, local_frame_xyz);

  int usage = 0;
  usage |= gpu::use_xyz;
  usage |= gpu::use_energy;
  usage |= gpu::use_grad;
  usage |= gpu::use_virial;

  SECTION("epolar -- pme") {
    std::string key1 = key0;
    key1 += "polarizeterm    only\n";
    key1 += "ewald\n";
    key1 += "ewald-cutoff    7.0\n";
    key1 += "neighbor-list\n";
    key1 += "list-buffer    0.1\n";
    key1 += "a-axis    20.0\n";
    file fke(k, key1);

    const char* argv[] = {"dummy", x1};
    int argc = 2;

    const double eps_e = 0.0001;
    // const double ref_eng = -36.5477;
    const double ref_eng = 0;
    const int ref_count = 222;

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    gpu::zero_egv();
    gpu::elec_init(gpu::v0);
    tinker_gpu_epolar0();
    gpu::torque(gpu::v0);
    COMPARE_ENERGY_(gpu::ep, ref_eng, eps_e);

    tinker_gpu_data_destroy();
    test_end();
  }
}
