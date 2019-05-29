#include "files.h"
#include "test/ff.h"
#include "test/rt.h"
#include "test/test.h"

m_tinker_using_namespace;
using namespace test;

#define COMPARE_CODE_BLOCK1_                                                   \
  {                                                                            \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal0();                                                    \
    COMPARE_ENERGY_(gpu::ev, ref_eng, eps);                                    \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal1();                                                    \
    COMPARE_ENERGY_(gpu::ev, ref_eng, eps);                                    \
    COMPARE_GRADIENT_(ref_grad, eps);                                          \
    COMPARE_VIR_(gpu::vir_ev, ref_v, eps);                                     \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal3();                                                    \
    COMPARE_ENERGY_(gpu::ev, ref_eng, eps);                                    \
    COMPARE_COUNT_(gpu::nev, ref_count);                                       \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal4();                                                    \
    COMPARE_ENERGY_(gpu::ev, ref_eng, eps);                                    \
    COMPARE_GRADIENT_(ref_grad, eps);                                          \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal5();                                                    \
    COMPARE_GRADIENT_(ref_grad, eps);                                          \
                                                                               \
    gpu::zero_egv();                                                           \
    tinker_gpu_evdw_hal6();                                                    \
    COMPARE_GRADIENT_(ref_grad, eps);                                          \
    COMPARE_VIR_(gpu::vir_ev, ref_v, eps);                                     \
  }

TEST_CASE("NaCl-1", "[ff][evdw][hal][switch][nacl]") {
  file fpr("amoeba09.prm", amoeba09_prm);

  std::string key = nacl_key;
  key += "vdwterm    only\n";
  file fke("test_nacl.key", key);

  int usage = 0;
  usage |= gpu::use_xyz;
  usage |= gpu::use_energy;
  usage |= gpu::use_grad;
  usage |= gpu::use_virial;

  const double eps = 1.0e-5;

  SECTION("ehal -- no switch") {
    const char* x1 = "test_nacl.xyz";
    file fx1(x1, nacl_xyz1);

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

  SECTION("ehal -- switch, near cut") {
    const char* x2 = "test_nacl.xyz_2";
    file fx2(x2, nacl_xyz2);
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

  SECTION("ehal -- switch, near off") {
    const char* x3 = "test_nacl.xyz_3";
    file fx3(x3, nacl_xyz3);

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

TEST_CASE("NaCl-2", "[ff][empole][coulomb][nacl]") {
  file fpr("amoeba09.prm", amoeba09_prm);

  std::string key = nacl_key;
  key += "multipoleterm    only\n";
  file fke("test_nacl.key", key);

  int usage = 0;
  usage |= gpu::use_xyz;
  usage |= gpu::use_energy;
  usage |= gpu::use_grad;
  usage |= gpu::use_virial;

  const double eps = 1.0e-5;

  SECTION("empole -- coulomb, pbc") {
    const char* x1 = "test_nacl.xyz";
    file fx1(x1, nacl_xyz1);

    const char* argv[] = {"dummy", x1};
    int argc = 2;

    const double ref_eng = -150.9381;
    const int ref_count = 1;
    const double ref_grad[][3] = {{-68.6082, 0.0, 0.0}, {68.6082, 0.0, 0.0}};
    const double ref_v[][3] = {
        {150.938, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    gpu::zero_egv();
    tinker_gpu_empole0();
    COMPARE_ENERGY_(gpu::em, ref_eng, eps);

    gpu::zero_egv();
    tinker_gpu_empole1();
    COMPARE_ENERGY_(gpu::em, ref_eng, eps);
    COMPARE_GRADIENT_(ref_grad, eps);
    COMPARE_VIR_(gpu::vir_em, ref_v, eps);

    gpu::zero_egv();
    tinker_gpu_empole3();
    COMPARE_ENERGY_(gpu::em, ref_eng, eps);
    COMPARE_COUNT_(gpu::nem, ref_count);

    tinker_gpu_data_destroy();
    test_end();
  }
}

TEST_CASE("NaCl-3", "[ff][empole][ewald][nacl]") {
  file fpr("amoeba09.prm", amoeba09_prm);

  std::string key = nacl_key4;
  key += "multipoleterm    only\n";
  file fke("test_nacl.key", key);

  int usage = 0;
  usage |= gpu::use_xyz;
  usage |= gpu::use_energy;
  usage |= gpu::use_grad;
  usage |= gpu::use_virial;

  const double eps = 1.0e-5;

  SECTION("empole -- pme") {
    const char* x4 = "test_nacl.xyz_4";
    file fx1(x4, nacl_xyz4);

    const char* argv[] = {"dummy", x4};
    int argc = 2;

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    gpu::zero_egv();
    tinker_gpu_empole0();

    tinker_gpu_data_destroy();
    test_end();
  }
}

#undef COMPARE_CODE_BLOCK1_
