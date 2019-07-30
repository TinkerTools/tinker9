#include "files.h"
#include "gpu/decl_md.h"
#include "test/ff.h"
#include "test/rt.h"
#include "test/test.h"

using namespace TINKER_NAMESPACE;
using namespace test;

static const char* verlet_intg = "integrator  verlet\n";
static int usage_ =
    gpu::use_xyz | gpu::use_vel | gpu::use_mass | gpu::vmask | gpu::use_md;

TEST_CASE("Kinetic-ArBox", "[ff][kinetic][arbox]") {
  const char* k = "test_arbox.key";
  const char* d = "test_arbox.dyn";
  const char* x = "test_arbox.xyz";
  const char* p = "amoeba09.prm";

  std::string k0 = arbox_key;
  k0 += verlet_intg;
  file fke(k, k0);

  file fd(d, arbox_dyn2);
  file fx(x, arbox_xyz);
  file fp(p, amoeba09_prm);

  const char* argv[] = {"dummy", x};
  int argc = 2;
  test_begin_1_xyz(argc, argv);
  test_mdinit(0, 0);

  use_data = usage_;
  tinker_gpu_runtime_initialize();

  real temp;
  gpu::kinetic(temp);

  const double ref_eksum = 100446.40376;
  const double ref_temp = 156008.001336;
  const double eps_e = 0.0001;

  REQUIRE(eksum == Approx(ref_eksum).margin(eps_e));
  REQUIRE(temp == Approx(ref_temp).margin(eps_e));

  tinker_gpu_runtime_finish();
  test_end();
}
