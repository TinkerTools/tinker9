#include "files.h"
#include "test/ff.h"
#include "test/rt.h"
#include "test/test.h"

m_tinker_using_namespace;
using namespace test;

static const char* triclinic_box = R"**(
ewald
ewald-cutoff  7.0
vdw-cutoff    7.0
neighbor-list
list-buffer  0.01

a-axis   28.0
b-axis   22.0
c-axis   18.0
alpha   105.0
beta    110.0
gamma    80.0
)**";

static int usage =
    gpu::use_xyz | gpu::use_energy | gpu::use_grad | gpu::use_virial;
static const char* k = "test_local_frame2.key";
static const char* x1 = "test_local_frame2.xyz";
static const char* argv[] = {"dummy", x1};
static int argc = 2;

TEST_CASE("Angle-Local-Frame2", "[ff][eangle][local-frame2]") {
  std::string k0 = local_frame_key;
  k0 += triclinic_box;
  k0 += "angleterm  only\n";
  file fke(k, k0);

  file fpr("amoeba09.prm", amoeba09_prm);
  file fx1(x1, local_frame_xyz2);

  test_begin_1_xyz(argc, argv);
  gpu::use_data = usage;
  tinker_gpu_data_create();

  const double eps_e = 0.0001;
  const double ref_e = 14.5726;
  const int ref_count = 12;
  const double eps_g = 0.0001;
  const double ref_g[][3] = {
      {0.0000, 0.0000, 0.0000},      {0.0000, 0.0000, 0.0000},
      {-5.7505, -1.8831, 1.4757},    {3.7651, -0.7057, 0.3221},
      {1.9854, 2.5888, -1.7977},     {-2.7022, 23.0671, -18.2046},
      {5.8412, -18.1715, -0.5945},   {-3.1390, -4.8955, 18.7991},
      {-7.1117, 1.0387, 8.1616},     {4.0056, -3.4254, -5.5430},
      {16.7852, -0.7352, 6.0062},    {-16.7783, 8.0355, -1.3414},
      {3.6707, -5.9216, -0.3892},    {-0.5714, 1.0081, -6.8942},
      {-9.0423, -12.1208, -53.3700}, {10.2396, 11.7722, 33.7171},
      {19.2501, -5.1824, 20.7327},   {-20.4475, 5.5310, -1.0798}};
  const double eps_v = 0.001;
  const double ref_v[][3] = {{-24.302, -5.249, -7.344},
                             {-5.249, 25.867, -23.021},
                             {-7.344, -23.021, -1.565}};

  COMPARE_BONED_FORCE(gpu::eangle, gpu::ea, ref_e, eps_e, gpu::nangle,
                      ref_count, gpu::gx, gpu::gy, gpu::gz, ref_g, eps_g,
                      gpu::vir_ea, ref_v, eps_v);

  tinker_gpu_data_destroy();
  test_end();
}
