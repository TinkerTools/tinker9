#include "files.h"
#include "test/ff.h"
#include "test/rt.h"
#include "test/test.h"

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

  file fx1(x1, nacl_xyz1);
  file fx2(x2, nacl_xyz2);
  file fx3(x3, nacl_xyz3);
  file fpr(prm, amoeba09_prm);
  file fke(k, key);

  int usage = 0;
  usage |= gpu::use_xyz;
  usage |= gpu::use_energy;
  usage |= gpu::use_grad;
  usage |= gpu::use_virial;

  const double eps = 1.0e-5;

  SECTION("no-switch") {
    const char* argv[] = {"dummy", x1};
    int argc = 2;

    const double ref_eng = 51.4242;
    const int ref_count = 1;
    const double ref_grad[][3] = {{184.4899, 0.0, 0.0}, {-184.4899, 0.0, 0.0}};
    const double ref_v[][3] = {
        {-405.878, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

#define COMPARE_ENERGY_                                                        \
  {                                                                            \
    double eng;                                                                \
    eng = gpu::get_energy(gpu::ev);                                            \
    REQUIRE(eng == Approx(ref_eng).epsilon(eps));                              \
  }
#define COMPARE_COUNT_                                                         \
  {                                                                            \
    int count;                                                                 \
    count = gpu::count_evdw();                                                 \
    REQUIRE(count == ref_count);                                               \
  }
#define COMPARE_GRAD_                                                          \
  {                                                                            \
    grad_t grad(gpu::n);                                                       \
    double* dst = &grad[0][0];                                                 \
    gpu::copyout_data2(0, 3, dst, gpu::gx, gpu::n);                            \
    gpu::copyout_data2(1, 3, dst, gpu::gy, gpu::n);                            \
    gpu::copyout_data2(2, 3, dst, gpu::gz, gpu::n);                            \
    for (int i = 0; i < gpu::n; ++i) {                                         \
      for (int j = 0; j < 3; ++j) {                                            \
        REQUIRE(grad[i][j] == Approx(ref_grad[i][j]).epsilon(eps));            \
      }                                                                        \
    }                                                                          \
  }
#define COMPARE_VIR_                                                           \
  {                                                                            \
    double vir[9];                                                             \
    gpu::get_virial(vir, gpu::vir_ev);                                         \
    for (int i = 0; i < 3; ++i) {                                              \
      for (int j = 0; j < 3; ++j) {                                            \
        int k = 3 * i + j;                                                     \
        REQUIRE(vir[k] == Approx(ref_v[i][j]).epsilon(eps));                   \
      }                                                                        \
    }                                                                          \
  }
#define COMPARE_CODE_BLOCK1_                                                   \
  {                                                                            \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal0();                                                    \
    COMPARE_ENERGY_;                                                           \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal1();                                                    \
    COMPARE_ENERGY_;                                                           \
    COMPARE_GRAD_;                                                             \
    COMPARE_VIR_;                                                              \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal3();                                                    \
    COMPARE_ENERGY_;                                                           \
    COMPARE_COUNT_;                                                            \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal4();                                                    \
    COMPARE_ENERGY_;                                                           \
    COMPARE_GRAD_;                                                             \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal5();                                                    \
    COMPARE_GRAD_;                                                             \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal6();                                                    \
    COMPARE_GRAD_;                                                             \
    COMPARE_VIR_;                                                              \
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
    const double ref_v[][3] = {
        {-354.835, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

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
    const double ref_v[][3] = {
        {-319.160, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_data_destroy();
    test_end();
  }
}

#undef COMPARE_ENERGY_
#undef COMPARE_COUNT_
#undef COMPARE_GRAD_
#undef COMPARE_VIR_
#undef COMPARE_CODE_BLOCK1_
