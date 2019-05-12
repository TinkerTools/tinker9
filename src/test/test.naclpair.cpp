#include "files.h"
#include "test/test.ff.h"
#include "test/test.h"
#include "test/test.os.h"

m_tinker_using_namespace;
using namespace test;

TEST_CASE("EHal-Switch-NaCl", "[forcefield]") {
  const char* x1 = "test_nacl.xyz";
  const char* x2 = "test_nacl.xyz_2";
  const char* x3 = "test_nacl.xyz_3";
  const char* prm = "amoeba09.prm";
  const char* k = "test_nacl.key";

  std::string key = nacl_key;
  key += "multipoleterm                  none\n";
  key += "polarizeterm                   none\n";

  file_gen fx1(x1, nacl_xyz1);
  file_gen fx2(x2, nacl_xyz2);
  file_gen fx3(x3, nacl_xyz3);
  file_gen fpr(prm, amoeba09_prm);
  file_gen fke(k, key);

  int usage = 0;
  usage |= use_xyz;
  usage |= use_energy;
  usage |= use_grad;
  usage |= use_virial;

  const double eps = 1.0e-5;

  SECTION("no-switch") {
    const char* argv[] = {"dummy", x1};
    int argc = 2;

    const double ref_eng = 51.4242;
    const int ref_count = 1;
    const double ref_grad[][3] = {{184.4899, 0.0, 0.0}, {-184.4899, 0.0, 0.0}};

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

#define COMPARE_CODE_BLOCK1_                                                   \
  {                                                                            \
    double eng;                                                                \
    int count;                                                                 \
    grad_t grad(gpu::n);                                                       \
    double* dst = &grad[0][0];                                                 \
                                                                               \
    tinker_gpu_zero_vag();                                                     \
    tinker_gpu_evdw_hal0();                                                    \
    eng = gpu::get_evdw();                                                     \
    REQUIRE(eng == Approx(ref_eng).epsilon(eps));                              \
                                                                               \
    tinker_gpu_zero_vag();                                                     \
    tinker_gpu_evdw_hal1();                                                    \
    eng = gpu::get_evdw();                                                     \
    REQUIRE(eng == Approx(ref_eng).epsilon(eps));                              \
    gpu::copyout_data_n(0, 3, dst, gpu::gx, gpu::n);                           \
    gpu::copyout_data_n(1, 3, dst, gpu::gy, gpu::n);                           \
    gpu::copyout_data_n(2, 3, dst, gpu::gz, gpu::n);                           \
    for (int i = 0; i < gpu::n; ++i)                                           \
      for (int j = 0; j < 3; ++j) {                                            \
        REQUIRE(grad[i][j] == Approx(ref_grad[i][j]).epsilon(eps));            \
      }                                                                        \
                                                                               \
    tinker_gpu_zero_vag();                                                     \
    tinker_gpu_evdw_hal3();                                                    \
    eng = gpu::get_evdw();                                                     \
    count = gpu::count_evdw();                                                 \
    REQUIRE(eng == Approx(ref_eng).epsilon(eps));                              \
    REQUIRE(count == ref_count);                                               \
  }

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_data_destroy();
    test_end();
  }

  SECTION("near-cut") {
    const char* argv[] = {"dummy", x2};
    int argc = 2;

    const double ref_eng = 25.8420;
    const int ref_count = 1;
    const double ref_grad[][3] = {{149.0904, 0.0, 0.0}, {-149.0904, 0.0, 0.0}};

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_data_destroy();
    test_end();
  }

  SECTION("near-off") {
    const char* argv[] = {"dummy", x3};
    int argc = 2;

    const double ref_eng = 4.8849;
    const int ref_count = 1;
    const double ref_grad[][3] = {{127.6639, 0.0, 0.0}, {-127.6639, 0.0, 0.0}};

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_data_destroy();
    test_end();
  }
}

#undef COMPARE_CODE_BLOCK1_
