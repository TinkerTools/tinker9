#include "files.h"
#include "test/ff.h"
#include "test/rt.h"
#include "test/test.h"

m_tinker_using_namespace;
using namespace test;

TEST_CASE("Local_Frame_1", "[forcefield][empole][coulomb][local_frame]") {
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

  const double eps = 1.0e-5;

  SECTION("case 1") {
    std::string key1 = key0;
    key1 += "multipoleterm    only\n";
    file fke(k, key1);

    const char* argv[] = {"dummy", x1};
    int argc = 2;

    const double ref_eng = -95.7432;
    const int ref_count = 201;

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    gpu::zero_egv();
    tinker_gpu_empole0();
    COMPARE_ENERGY_(gpu::em, ref_eng, eps);

    gpu::zero_egv();
    tinker_gpu_empole3();
    COMPARE_ENERGY_(gpu::em, ref_eng, eps);
    COMPARE_COUNT_(gpu::nem, ref_count);

    tinker_gpu_data_destroy();
    test_end();
  }
}
