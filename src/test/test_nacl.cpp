#include "files.h"
#include "test/ff.h"
#include "test/rt.h"
#include "test/test.h"

m_tinker_using_namespace;
using namespace test;

#define COMPARE_ENERGY_(gpuptr)                                                \
  {                                                                            \
    double eng = gpu::get_energy(gpuptr);                                      \
    REQUIRE(eng == Approx(ref_eng).epsilon(eps));                              \
  }
#define COMPARE_COUNT_(gpuptr)                                                 \
  {                                                                            \
    int count = gpu::get_count(gpuptr);                                        \
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
#define COMPARE_VIR_(gpuptr)                                                   \
  {                                                                            \
    double vir[9];                                                             \
    gpu::get_virial(vir, gpuptr);                                              \
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
    COMPARE_ENERGY_(gpu::ev);                                                  \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal1();                                                    \
    COMPARE_ENERGY_(gpu::ev);                                                  \
    COMPARE_GRAD_;                                                             \
    COMPARE_VIR_(gpu::vir_ev);                                                 \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal3();                                                    \
    COMPARE_ENERGY_(gpu::ev);                                                  \
    COMPARE_COUNT_(gpu::nev);                                                  \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal4();                                                    \
    COMPARE_ENERGY_(gpu::ev);                                                  \
    COMPARE_GRAD_;                                                             \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal5();                                                    \
    COMPARE_GRAD_;                                                             \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal6();                                                    \
    COMPARE_GRAD_;                                                             \
    COMPARE_VIR_(gpu::vir_ev);                                                 \
  }

TEST_CASE("EHal-Switch-NaCl", "[forcefield][ehal][nacl]") {
  const char* x1 = "test_nacl.xyz";
  const char* x2 = "test_nacl.xyz_2";
  const char* x3 = "test_nacl.xyz_3";
  const char* prm = "amoeba09.prm";
  const char* k = "test_nacl.key";

  std::string key = nacl_key;
  key += "vdwterm                        only\n";

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

  SECTION("case 1, no-switch") {
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

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_data_destroy();
    test_end();
  }

  SECTION("case 2, near-cut") {
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

  SECTION("case 3, near-off") {
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

TEST_CASE("EMpole-Coulomb-NaCl", "[forcefield][empole][nacl]") {
  const char* x1 = "test_nacl.xyz";
  const char* prm = "amoeba09.prm";
  const char* k = "test_nacl.key";

  std::string key = nacl_key;
  key += "multipoleterm                  only\n";

  file fx1(x1, nacl_xyz1);
  file fpr(prm, amoeba09_prm);
  file fke(k, key);

  int usage = 0;
  usage |= gpu::use_xyz;
  usage |= gpu::use_energy;
  usage |= gpu::use_grad;
  usage |= gpu::use_virial;

  const double eps = 1.0e-5;

  SECTION("case 1") {
    const char* argv[] = {"dummy", x1};
    int argc = 2;

    const double ref_eng = -150.9381;
    const int ref_count = 1;

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    gpu::zero_egv();
    tinker_gpu_empole0();
    COMPARE_ENERGY_(gpu::em);

    gpu::zero_egv();
    tinker_gpu_empole3();
    COMPARE_ENERGY_(gpu::em);
    COMPARE_COUNT_(gpu::nem);

    tinker_gpu_data_destroy();
    test_end();
  }
}

#undef COMPARE_ENERGY_
#undef COMPARE_COUNT_
#undef COMPARE_GRAD_
#undef COMPARE_VIR_
#undef COMPARE_CODE_BLOCK1_
