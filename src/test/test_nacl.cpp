#include "files.h"
#include "test/ff.h"
#include "test/rt.h"
#include "test/test.h"

using namespace TINKER_NAMESPACE;
using namespace test;

#define COMPARE_CODE_BLOCK1_                                                   \
  {                                                                            \
    gpu::zero_egv();                                                           \
    gpu::evdw_hal(gpu::v0);                                                    \
    COMPARE_ENERGY_(gpu::ev, ref_eng, eps);                                    \
                                                                               \
    gpu::zero_egv();                                                           \
    gpu::evdw_hal(gpu::v1);                                                    \
    COMPARE_ENERGY_(gpu::ev, ref_eng, eps);                                    \
    COMPARE_GRADIENT_(ref_grad, eps);                                          \
    COMPARE_VIR_(gpu::vir_ev, ref_v, eps);                                     \
                                                                               \
    gpu::zero_egv();                                                           \
    gpu::evdw_hal(gpu::v3);                                                    \
    COMPARE_ENERGY_(gpu::ev, ref_eng, eps);                                    \
    COMPARE_COUNT_(gpu::nev, ref_count);                                       \
                                                                               \
    gpu::zero_egv();                                                           \
    gpu::evdw_hal(gpu::v4);                                                    \
    COMPARE_ENERGY_(gpu::ev, ref_eng, eps);                                    \
    COMPARE_GRADIENT_(ref_grad, eps);                                          \
                                                                               \
    gpu::zero_egv();                                                           \
    gpu::evdw_hal(gpu::v5);                                                    \
    COMPARE_GRADIENT_(ref_grad, eps);                                          \
                                                                               \
    gpu::zero_egv();                                                           \
    gpu::evdw_hal(gpu::v6);                                                    \
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
  usage |= gpu::vmask;

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
    tinker_gpu_runtime_initialize();

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_runtime_finish();
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
    tinker_gpu_runtime_initialize();

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_runtime_finish();
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
    tinker_gpu_runtime_initialize();

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_runtime_finish();
    test_end();
  }
}

#define COMPARE_CODE_BLOCK2_                                                   \
  {                                                                            \
    gpu::zero_egv();                                                           \
    gpu::elec_init(gpu::v0);                                                   \
    gpu::empole(gpu::v0);                                                      \
    gpu::torque(gpu::v0);                                                      \
    COMPARE_ENERGY_(gpu::em, ref_eng, eps);                                    \
                                                                               \
    gpu::zero_egv();                                                           \
    gpu::elec_init(gpu::v1);                                                   \
    gpu::empole(gpu::v1);                                                      \
    gpu::torque(gpu::v1);                                                      \
    COMPARE_ENERGY_(gpu::em, ref_eng, eps);                                    \
    COMPARE_GRADIENT_(ref_grad, eps_g);                                        \
    COMPARE_VIR2_(gpu::vir_em, gpu::vir_trq, ref_v, eps_v);                    \
                                                                               \
    gpu::zero_egv();                                                           \
    gpu::elec_init(gpu::v3);                                                   \
    gpu::empole(gpu::v3);                                                      \
    gpu::torque(gpu::v3);                                                      \
    COMPARE_ENERGY_(gpu::em, ref_eng, eps);                                    \
    COMPARE_COUNT_(gpu::nem, ref_count);                                       \
                                                                               \
    gpu::zero_egv();                                                           \
    gpu::elec_init(gpu::v4);                                                   \
    gpu::empole(gpu::v4);                                                      \
    gpu::torque(gpu::v4);                                                      \
    COMPARE_ENERGY_(gpu::em, ref_eng, eps);                                    \
    COMPARE_GRADIENT_(ref_grad, eps_g);                                        \
                                                                               \
    gpu::zero_egv();                                                           \
    gpu::elec_init(gpu::v5);                                                   \
    gpu::empole(gpu::v5);                                                      \
    gpu::torque(gpu::v5);                                                      \
    COMPARE_GRADIENT_(ref_grad, eps_g);                                        \
                                                                               \
    gpu::zero_egv();                                                           \
    gpu::elec_init(gpu::v6);                                                   \
    gpu::empole(gpu::v6);                                                      \
    gpu::torque(gpu::v6);                                                      \
    COMPARE_GRADIENT_(ref_grad, eps_g);                                        \
    COMPARE_VIR2_(gpu::vir_em, gpu::vir_trq, ref_v, eps_v);                    \
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
    const double eps_g = eps;
    const double ref_grad[][3] = {{-68.6082, 0.0, 0.0}, {68.6082, 0.0, 0.0}};
    const double eps_v = eps;
    const double ref_v[][3] = {
        {150.938, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_runtime_initialize();

    COMPARE_CODE_BLOCK2_;

    tinker_gpu_runtime_finish();
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

  const double eps = 5.0e-4;

  SECTION("empole -- pme") {
    const char* x4 = "test_nacl.xyz_4";
    file fx1(x4, nacl_xyz4);

    const char* argv[] = {"dummy", x4};
    int argc = 2;

    const double ref_eng = -0.3146;
    const int ref_count = 2;
    // total grad
    const double eps_g = eps;
    const double ref_grad[][3] = {{0.1731, 0.1921, 0.2103},
                                  {-0.1668, -0.1917, -0.2081}};
    // self grad = 0
    // real grad
    // const double ref_grad[][3] = {{21.9459, 24.1404, 26.3789},
    //                               {-21.9459, -24.1404, -26.3789}};
    // recip grad
    // const double ref_grad[][3] = {{-21.7728, -23.9484, -26.1686},
    //                               {21.7791, 23.9487, 26.1709}};
    const double eps_v = eps;
    // total virial
    const double ref_v[][3] = {{0.059, -0.284, -0.311},
                               {-0.284, 0.105, -0.342},
                               {-0.311, -0.342, 0.155}};
    // self virial = 0
    // real virial
    // const double ref_v[][3] = {{-21.946, -24.140, -26.379},
    //                            {-24.140, -26.554, -29.017},
    //                            {-26.379, -29.017, -31.707}};
    // recip virial
    // const double ref_v[][3] = {{22.005, 23.857, 26.068},
    //                            {23.857, 26.659, 28.675},
    //                            {26.068, 28.675, 31.863}};

    test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_runtime_initialize();

    COMPARE_CODE_BLOCK2_;

    tinker_gpu_runtime_finish();
    test_end();
  }
}

#undef COMPARE_CODE_BLOCK1_
#undef COMPARE_CODE_BLOCK2_
