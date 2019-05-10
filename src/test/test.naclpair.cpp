#include "files.h"
#include "gpu/data.h"
#include "gpu/mdstate.h"
#include "gpu/potential.h"
#include "test/test.h"
#include "test/test.os.h"

m_tinker_using_namespace;

TEST_CASE("EHal-Switch-NaCl", "[forcefield]") {
  const char* x1 = "test_nacl.xyz";
  const char* x2 = "test_nacl.xyz_2";
  const char* x3 = "test_nacl.xyz_3";
  const char* prm = "amoeba09.prm";
  const char* k = "test_nacl.key";

  std::string key = test::nacl_key;
  key += "multipoleterm                  none\n";
  key += "polarizeterm                   none\n";

  test::file_gen fx1(x1, test::nacl_xyz1);
  test::file_gen fx2(x2, test::nacl_xyz2);
  test::file_gen fx3(x3, test::nacl_xyz3);
  test::file_gen fpr(prm, test::amoeba09_prm);
  test::file_gen fke(k, key);

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

    test::test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

#define COMPARE_CODE_BLOCK1_                                                   \
  {                                                                            \
    double eng;                                                                \
    int count;                                                                 \
    tinker_gpu_evdw_hal1();                                                    \
    eng = gpu::get_evdw();                                                     \
    REQUIRE(eng == Approx(ref_eng).epsilon(eps));                              \
                                                                               \
    tinker_gpu_evdw_hal3();                                                    \
    eng = gpu::get_evdw();                                                     \
    count = gpu::count_evdw();                                                 \
    REQUIRE(eng == Approx(ref_eng).epsilon(eps));                              \
    REQUIRE(count == ref_count);                                               \
  }

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_data_destroy();
    test::test_end();
  }

  SECTION("near-cut") {
    const char* argv[] = {"dummy", x2};
    int argc = 2;

    const double ref_eng = 25.8420;
    const int ref_count = 1;

    test::test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_data_destroy();
    test::test_end();
  }

  SECTION("near-off") {
    const char* argv[] = {"dummy", x3};
    int argc = 2;

    const double ref_eng = 4.8849;
    const int ref_count = 1;

    test::test_begin_1_xyz(argc, argv);
    gpu::use_data = usage;
    tinker_gpu_data_create();

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_data_destroy();
    test::test_end();
  }
}

#undef COMPARE_CODE_BLOCK1_
