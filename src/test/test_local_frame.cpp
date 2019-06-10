#include "files.h"
#include "test/ff.h"
#include "test/rt.h"
#include "test/test.h"
#include <ext/tinker/tinker_mod.h>

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
    key1 += "usolve-cutoff    0.01\n";
    key1 += "polarizeterm    only\n";
    key1 += "ewald\n";
    key1 += "ewald-cutoff    7.0\n";
    key1 += "neighbor-list\n";
    key1 += "list-buffer    0.1\n";
    key1 += "a-axis    20.0\n";
    file fke(k, key1);

    const char* argv[] = {"dummy", x1};
    int argc = 2;

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    // dfield ewald version

    const double eps_f = 0.0001;
    const double ref_dir_field_d[][3] = {
        {0.0603, -0.0877, -0.0024}, {0.0304, -0.0254, 0.0053},
        {0.0743, -0.0256, -0.0002}, {0.0540, -0.0080, -0.0020},
        {0.0577, -0.0410, 0.0126},  {0.1771, -0.1487, 0.0468},
        {0.1263, -0.1083, 0.0180},  {0.2235, -0.1782, 0.0293},
        {0.0589, -0.0228, 0.0400},  {0.0394, -0.0297, -0.0007},
        {0.0204, -0.0223, -0.0171}, {0.0144, -0.0136, -0.0118},
        {0.0160, -0.0218, -0.0204}, {0.0188, -0.0281, -0.0310},
        {-0.0052, -0.0557, 0.0866}, {0.0372, -0.0710, 0.0609},
        {-0.0302, -0.0249, 0.0959}, {-0.0002, -0.0971, 0.1573},
        {0.0176, -0.1618, -0.0356}, {0.0020, -0.0699, -0.0279},
        {0.1436, -0.1515, 0.0008},  {-0.1249, -0.1932, -0.0900}};
    const double ref_dir_field_p[][3] = {
        {0.0603, -0.0877, -0.0024}, {0.0304, -0.0254, 0.0053},
        {0.0743, -0.0256, -0.0002}, {0.0540, -0.0080, -0.0020},
        {0.0577, -0.0410, 0.0126},  {0.1771, -0.1487, 0.0468},
        {0.1263, -0.1083, 0.0180},  {0.2235, -0.1782, 0.0293},
        {0.0670, -0.0261, 0.0454},  {0.0394, -0.0297, -0.0007},
        {0.0204, -0.0223, -0.0171}, {0.0345, -0.0218, 0.0190},
        {0.0577, -0.0249, -0.0034}, {0.0456, -0.0554, -0.0137},
        {-0.0052, -0.0557, 0.0866}, {0.0372, -0.0710, 0.0609},
        {-0.0302, -0.0249, 0.0959}, {-0.0002, -0.0971, 0.1573},
        {0.0176, -0.1618, -0.0356}, {0.0020, -0.0699, -0.0279},
        {0.1436, -0.1515, 0.0008},  {-0.1249, -0.1932, -0.0900}};
    gpu::zero_egv();
    gpu::elec_init(gpu::v0);
    gpu::dfield_ewald(&gpu::udir[0][0], &gpu::udirp[0][0]);
    grad_t fieldd, fieldp;
    gpu::copyout_data3(fieldd, gpu::udir, gpu::n);
    gpu::copyout_data3(fieldp, gpu::udirp, gpu::n);
    for (int i = 0; i < gpu::n; ++i) {
      for (int j = 0; j < 3; ++j) {
        REQUIRE(fieldd[i][j] == Approx(ref_dir_field_d[i][j]).margin(eps_f));
        REQUIRE(fieldp[i][j] == Approx(ref_dir_field_p[i][j]).margin(eps_f));
      }
    }

    // ufield ewald version

    const double eps_uf = 0.0001;
    const double ref_ufield_d[][3] = {
        {0.2326, 1.7537, 0.4066},    {-0.0625, -0.7063, 0.1748},
        {-0.4040, -0.8194, -0.7421}, {-0.3634, -0.4587, -0.5382},
        {-0.3778, -1.1806, -0.6153}, {-0.3945, -0.4370, 0.1095},
        {-0.5747, -0.2189, 0.5939},  {-0.2782, -0.4352, 0.2659},
        {0.5975, -0.5023, -0.4645},  {0.1915, -1.2181, 0.0277},
        {-1.3381, -1.3854, -0.6932}, {-1.0270, -0.5388, 0.7093},
        {-0.0693, 0.0332, -0.1906},  {-1.0449, 0.0542, -0.1094},
        {-1.7474, -1.6330, -1.6222}, {-0.7535, -1.2197, -1.0019},
        {-1.2575, 0.2898, -0.0987},  {-0.9060, -1.1180, 0.6608},
        {-1.9284, -2.2206, -2.4186}, {0.0085, 0.3752, -1.6604},
        {0.7407, -1.4251, -1.0051},  {-1.4916, -0.8335, -1.4794}};
    const double ref_ufield_p[][3] = {
        {0.2088, 1.6503, 0.4078},    {-0.0632, -0.6470, 0.1741},
        {-0.3752, -0.6846, -0.5516}, {-0.3275, -0.3521, -0.3940},
        {-0.3701, -1.0143, -0.4764}, {-0.3450, -0.3801, 0.1822},
        {-0.4963, -0.2093, 0.5370},  {-0.2139, -0.4325, 0.2956},
        {0.5692, -0.4457, -0.4214},  {0.1598, -1.1058, 0.0538},
        {-1.2888, -1.2553, -0.5901}, {-0.9675, -0.4872, 0.6128},
        {-0.0813, 0.0744, -0.1494},  {-0.9896, 0.0318, -0.0839},
        {-1.6770, -1.5217, -1.4601}, {-0.7271, -1.1385, -0.8880},
        {-1.1935, 0.2491, -0.0607},  {-0.8772, -1.0243, 0.6095},
        {-1.8666, -2.1138, -2.2301}, {-0.0191, 0.4125, -1.5382},
        {0.7158, -1.3417, -0.9288},  {-1.4493, -0.8422, -1.3602}};
    gpu::zero_egv();
    gpu::elec_init(gpu::v0);
    grad_t ud, up;
    ud.resize(gpu::n);
    up.resize(gpu::n);
    for (int i = 0; i < gpu::n; ++i) {
      for (int j = 0; j < 3; ++j) {
        ud[i][j] = 0.1 * (i + 1) + 0.03 * (j + 1);
        up[i][j] = 0.1 * (i + 1) - 0.03 * (j + 1);
      }
    }
    gpu::copyin_data(&gpu::uind[0][0], &ud[0][0], 3 * gpu::n);
    gpu::copyin_data(&gpu::uinp[0][0], &up[0][0], 3 * gpu::n);
    gpu::ufield_ewald(&gpu::uind[0][0], &gpu::uinp[0][0], &gpu::udir[0][0],
                      &gpu::udirp[0][0]);
    gpu::copyout_data3(fieldd, gpu::udir, gpu::n);
    gpu::copyout_data3(fieldp, gpu::udirp, gpu::n);
    const double debye = units::debye;
    for (int i = 0; i < gpu::n; ++i) {
      for (int j = 0; j < 3; ++j) {
        REQUIRE(fieldd[i][j] == Approx(ref_ufield_d[i][j]).margin(eps_uf));
        REQUIRE(fieldp[i][j] == Approx(ref_ufield_p[i][j]).margin(eps_uf));
      }
    }

    // diagnoal matrix pcg

    const double eps_ind = 0.0001;
    const double ref_ud_debye[][3] = {
        {0.0243, -0.0845, -0.0109}, {1.2653, -0.8136, 0.5372},
        {0.2268, -0.0683, 0.0078},  {0.0616, -0.0298, 0.0188},
        {0.1019, -0.0927, 0.0913},  {0.5383, -0.5306, 0.1597},
        {0.1600, -0.1537, -0.0985}, {0.5290, -0.4213, -0.0113},
        {0.2702, -0.0662, 0.1661},  {0.2420, -0.1248, -0.0516},
        {0.1211, -0.0628, -0.0252}, {0.0168, -0.0089, -0.0295},
        {0.0499, -0.0125, -0.0515}, {0.0365, -0.0617, -0.0396},
        {-0.0212, -0.1500, 0.3140}, {0.1487, -0.1622, 0.1192},
        {-0.0466, 0.0288, 0.1480},  {-0.0021, -0.1708, 0.4490},
        {0.2273, -0.6081, -0.0940}, {0.0143, -0.1151, -0.0139},
        {0.3569, -0.1912, 0.0563},  {-0.2176, -0.4640, -0.1812}};
    const double ref_up_debye[][3] = {
        {0.0241, -0.0843, -0.0110}, {1.2706, -0.8115, 0.5494},
        {0.2253, -0.0684, 0.0097},  {0.0612, -0.0301, 0.0192},
        {0.1011, -0.0925, 0.0915},  {0.5377, -0.5296, 0.1603},
        {0.1596, -0.1533, -0.0979}, {0.5288, -0.4205, -0.0101},
        {0.3149, -0.0836, 0.1982},  {0.2894, -0.1421, -0.0248},
        {0.0447, -0.0298, -0.0822}, {0.0657, -0.0302, 0.0391},
        {0.1486, -0.0288, -0.0009}, {0.1037, -0.1202, 0.0087},
        {-0.0179, -0.1459, 0.3136}, {0.1491, -0.1603, 0.1190},
        {-0.0422, 0.0305, 0.1471},  {-0.0008, -0.1699, 0.4506},
        {0.2227, -0.6054, -0.0918}, {0.0119, -0.1147, -0.0138},
        {0.3551, -0.1907, 0.0577},  {-0.2191, -0.4626, -0.1832}};
    gpu::zero_egv();
    gpu::elec_init(gpu::v0);
    gpu::induce(&gpu::uind[0][0], &gpu::uinp[0][0]);
    gpu::copyout_data3(ud, gpu::uind, gpu::n);
    gpu::copyout_data3(up, gpu::uinp, gpu::n);
    for (int i = 0; i < gpu::n; ++i) {
      for (int j = 0; j < 3; ++j) {
        REQUIRE(ud[i][j] * debye == Approx(ref_ud_debye[i][j]).margin(eps_ind));
        REQUIRE(up[i][j] * debye == Approx(ref_up_debye[i][j]).margin(eps_ind));
      }
    }

    const double eps_e = 0.0001;
    // const double ref_eng = -36.5477;
    const double ref_eng = 0;
    const int ref_count = 222;

    gpu::zero_egv();
    gpu::elec_init(gpu::v0);
    tinker_gpu_epolar0();
    gpu::torque(gpu::v0);
    COMPARE_ENERGY_(gpu::ep, ref_eng, eps_e);

    tinker_gpu_data_destroy();
    test_end();
  }
}
