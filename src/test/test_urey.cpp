#include "files.h"
#include "test.h"
#include "test_rt.h"

using namespace TINKER_NAMESPACE;

static const char* ureyterm_only = R"**(
ureyterm  only
)**";

static int usage = calc::xyz | calc::vmask;

static const double ref_g_urey_h2o10[][3] = {
    {0.0000, 0.0000, 0.0000},   {-0.1502, -0.0096, 0.2033},
    {0.1502, 0.0096, -0.2033},  {0.0000, 0.0000, 0.0000},
    {-0.1258, 0.0384, -0.8188}, {0.1258, -0.0384, 0.8188},
    {0.0000, 0.0000, 0.0000},   {0.0444, -0.0151, 0.0147},
    {-0.0444, 0.0151, -0.0147}, {0.0000, 0.0000, 0.0000},
    {-0.2235, 1.0113, 0.4196},  {0.2235, -1.0113, -0.4196},
    {0.0000, 0.0000, 0.0000},   {-0.3173, -0.0814, 0.6494},
    {0.3173, 0.0814, -0.6494},  {0.0000, 0.0000, 0.0000},
    {0.6837, -0.4755, -0.4068}, {-0.6837, 0.4755, 0.4068},
    {0.0000, 0.0000, 0.0000},   {-0.2546, -0.7133, 0.1591},
    {0.2546, 0.7133, -0.1591},  {0.0000, 0.0000, 0.0000},
    {-0.2182, -0.0213, 0.1014}, {0.2182, 0.0213, -0.1014},
    {0.0000, 0.0000, 0.0000},   {-0.5657, 0.4544, -0.6075},
    {0.5657, -0.4544, 0.6075},  {0.0000, 0.0000, 0.0000},
    {0.2971, -0.0299, -0.9680}, {-0.2971, 0.0299, 0.9680},
};
TEST_CASE("Urey-Ten-Water", "[ff][eurey][h2o10]") {
  const char* k = "test_h2o10.key";
  const char* x1 = "test_h2o10.xyz";
  const char* p = "water03.prm";

  std::string k0 = h2o10_key;
  k0 += ureyterm_only;
  TestFile fke(k, k0);

  TestFile fx1(x1, h2o10_xyz);
  TestFile fpr(p, water03_prm);

  const char* argv[] = {"dummy", x1};
  int argc = 2;
  test_begin_with_args(argc, argv);
  rc_flag = usage;
  initialize();

  const double eps_e = 0.0001;
  const double ref_e = -0.1964;
  const int ref_count = 10;
  const double eps_g = 0.0001;
  const double eps_v = 0.001;
  const double ref_v[][3] = {
      {2.318, -0.827, -1.078}, {-0.827, 3.061, 0.062}, {-1.078, 0.062, 4.915}};

  COMPARE_BONDED_FORCE(eurey, eub, ref_e, eps_e, nurey, ref_count, gx, gy, gz,
                       ref_g_urey_h2o10, eps_g, vir_eub, ref_v, eps_v);

  finish();
  test_end();
}
