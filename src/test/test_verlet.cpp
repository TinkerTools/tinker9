#include "files.h"
#include "gpu/decl_md.h"
#include "test/ff.h"
#include "test/rt.h"
#include "test/test.h"

m_tinker_using_namespace;
using namespace test;

const char* verlet_intg = "integrator  verlet\n";
static int usage_ = gpu::use_xyz | gpu::use_vel | gpu::use_mass |
    gpu::use_energy | gpu::use_grad;

TEST_CASE("Verlet-Ar6", "[noassert][verlet][ar6]") {
  const char* k = "test_ar6.key";
  const char* d = "test_ar6.dyn";
  const char* x1 = "test_ar6.xyz";
  const char* p = "amoeba09.prm";

  std::string k0 = ar6_key;
  k0 += verlet_intg;
  file fke(k, k0);

  file fd(d, ar6_dyn);
  file fx1(x1, ar6_xyz);
  file fpr(p, amoeba09_prm);

  const char* argv[] = {"dummy", x1};
  int argc = 2;
  test_begin_1_xyz(argc, argv);
  test_mdinit(0, 0);

  gpu::use_data = usage_;
  tinker_gpu_data_create();

  gpu::md_data(static_cast<gpu::rc_t>(gpu::rc_alloc | gpu::rc_copyin));
  gpu::velocity_verlet(1, 0.001);

  tinker_gpu_data_destroy();
  test_end();
}
