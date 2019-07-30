#include "files.h"
#include "test/ff.h"
#include "test/rt.h"
#include "test/test.h"

using namespace TINKER_NAMESPACE;
using namespace test;

#define COMPARE_GRAD_                                                          \
  {                                                                            \
    grad_t grad(n);                                                            \
    double* dst = &grad[0][0];                                                 \
    copyout_array2(0, 3, dst, gx, n);                                          \
    copyout_array2(1, 3, dst, gy, n);                                          \
    copyout_array2(2, 3, dst, gz, n);                                          \
    for (int i = 0; i < 6; ++i) {                                              \
      for (int j = 0; j < 3; ++j) {                                            \
        REQUIRE(grad[i * 30][j] == Approx(ref_grad[i][j]).margin(eps));        \
      }                                                                        \
    }                                                                          \
  }
#define COMPARE_CODE_BLOCK1_                                                   \
  {                                                                            \
    zero_egv();                                                                \
    evdw_hal(calc::v0);                                                        \
    COMPARE_ENERGY_(ev, ref_eng, eps);                                         \
                                                                               \
    zero_egv();                                                                \
    evdw_hal(calc::v1);                                                        \
    COMPARE_ENERGY_(ev, ref_eng, eps);                                         \
    COMPARE_GRAD_;                                                             \
    COMPARE_VIR_(vir_ev, ref_v, eps_v);                                        \
                                                                               \
    zero_egv();                                                                \
    evdw_hal(calc::v3);                                                        \
    COMPARE_ENERGY_(ev, ref_eng, eps);                                         \
    COMPARE_COUNT_(nev, ref_count);                                            \
                                                                               \
    zero_egv();                                                                \
    evdw_hal(calc::v4);                                                        \
    COMPARE_ENERGY_(ev, ref_eng, eps);                                         \
    COMPARE_GRAD_;                                                             \
                                                                               \
    zero_egv();                                                                \
    evdw_hal(calc::v5);                                                        \
    COMPARE_GRAD_;                                                             \
                                                                               \
    zero_egv();                                                                \
    evdw_hal(calc::v6);                                                        \
    COMPARE_GRAD_;                                                             \
    COMPARE_VIR_(vir_ev, ref_v, eps_v);                                        \
  }

TEST_CASE("CLN025", "[ff][evdw][hal][cln025]") {
  const char* x = "test_cln025.xyz";
  const char* k = "test_cln025.key";
  const char* p = "amoebabio09.prm";
  std::string k0 = cln025_key;
  k0 += "vdwterm    only\n";

  file fx(x, cln025_xyz);
  file px(p, amoebabio09_prm);

  int usage = 0;
  usage |= calc::xyz;
  usage |= calc::vmask;

  SECTION("ehal -- gas phase, no cutoff") {
    file kx(k, k0);

    const char* argv[] = {"dummy", x};
    int argc = 2;

    const double eps = 2.0e-4;
    const double ref_eng = 88.8860;
    const int ref_count = 13225;
    // atom 1, 31, 61, 91, 121, 151
    const double ref_grad[][3] = {
        {-10.4202, 0.6937, 2.7019}, {4.0910, 2.0984, -2.4349},
        {-0.5266, 3.5665, 0.4037},  {-0.2577, -1.9100, 6.9858},
        {1.0470, 3.2855, -4.3019},  {4.3475, -2.2594, 1.8280}};
    const double eps_v = 0.001;
    const double ref_v[][3] = {{-800.488, -37.589, -2.250},
                               {-37.589, -758.657, 41.895},
                               {-2.250, 41.895, -681.179}};

    test_begin_1_xyz(argc, argv);
    use_data = usage;
    tinker_gpu_runtime_initialize();

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_runtime_finish();
    test_end();
  }

  SECTION("ehal -- pbc, cutoff") {
    std::string k1 = k0;
    k1 += "neighbor-list\n";
    k1 += "vdw-cutoff        6.0\n";
    k1 += "a-axis           18.0\n";
    file kx(k, k1);

    const char* argv[] = {"dummy", x};
    int argc = 2;

    const double eps = 2.0e-4;
    const double ref_eng = 104.3872;
    const int ref_count = 4116;
    // atom 1, 31, 61, 91, 121, 151
    const double ref_grad[][3] = {
        {-10.5674, 0.6411, 2.6721}, {4.0652, 2.1284, -2.4015},
        {-0.5214, 3.5720, 0.3980},  {-0.1419, -1.9235, 6.9387},
        {1.0707, 3.2989, -4.2483},  {4.3235, -2.2325, 1.8485}};
    const double eps_v = 0.001;
    const double ref_v[][3] = {{-914.956, -8.336, 26.361},
                               {-8.336, -781.555, 29.820},
                               {26.361, 29.820, -706.064}};

    test_begin_1_xyz(argc, argv);
    use_data = usage;
    tinker_gpu_runtime_initialize();

    COMPARE_CODE_BLOCK1_;

    tinker_gpu_runtime_finish();
    test_end();
  }
}

#undef COMPARE_GRAD_
#undef COMPARE_CODE_BLOCK1_
