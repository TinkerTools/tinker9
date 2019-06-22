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

TEST_CASE("Bond-Local-Frame2", "[ff][ebond][local-frame2]") {
  std::string k0 = local_frame_key;
  k0 += triclinic_box;
  k0 += "bondterm  only\n";
  file fke(k, k0);

  file fpr("amoeba09.prm", amoeba09_prm);
  file fx1(x1, local_frame_xyz2);

  test_begin_1_xyz(argc, argv);
  gpu::use_data = usage;
  tinker_gpu_data_create();

  const double eps_e = 0.0001;
  const double ref_e = 6.7067;
  const int ref_count = 12;
  const double eps_g = test_get_eps2(0.0002, 0.0001);
  const double ref_g[][3] = {
      {0.0000, 0.0000, 0.0000},      {0.0000, 0.0000, 0.0000},
      {9.0683, -22.1136, 14.3407},   {3.7174, 15.1673, -10.2237},
      {-12.7856, 6.9463, -4.1170},   {2.0528, 27.8731, -43.3919},
      {-10.3884, -4.9091, 47.9802},  {8.3356, -22.9641, -4.5883},
      {25.3020, -12.2233, 23.6028},  {-37.9336, 17.9922, -33.3279},
      {3.8409, -5.6886, 0.8513},     {-0.3951, 1.3165, 12.8289},
      {8.1264, 5.2279, -2.8986},     {1.0593, -6.6247, -1.0565},
      {38.2559, -79.9730, -52.4603}, {-36.9462, 38.1189, -2.0888},
      {1.8433, 44.7080, 9.4640},     {-3.1529, -2.8540, 45.0851}};
  const double eps_v = 0.001;
  const double ref_v[][3] = {{42.567, -52.373, 3.424},
                             {-52.373, -8.084, -13.810},
                             {3.424, -13.810, 85.120}};

  COMPARE_BONED_FORCE(gpu::ebond, gpu::eb, ref_e, eps_e, gpu::nbond, ref_count,
                      gpu::gx, gpu::gy, gpu::gz, ref_g, eps_g, gpu::vir_eb,
                      ref_v, eps_v);

  tinker_gpu_data_destroy();
  test_end();
}
