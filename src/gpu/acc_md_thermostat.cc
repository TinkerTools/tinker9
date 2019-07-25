#include "gpu/decl_mdstate.h"
#include "gpu/f_random.h"
#include <ext/tinker/tinker_mod.h>

TINKER_NAMESPACE_BEGIN
namespace gpu {
void thermo_bussi_acc_impl_(real dt_real, real temp_real) {
  double dt = dt_real;
  double temp = temp_real;

  double tautemp = bath::tautemp;
  double kelvin = bath::kelvin;
  int nfree = mdstuf::nfree;
  double& eta = bath::eta;

  if (temp == 0)
    temp = 0.1;

  double c = std::exp(-dt / tautemp);
  double d = (1 - c) * (kelvin / temp) / nfree;
  double r = normal_double();
  double s = chi_squared_double(nfree - 1);
  double scale = c + (s + r * r) * d + 2 * r * std::sqrt(c * d);
  scale = std::sqrt(scale);
  if (r + std::sqrt(c / d) < 0)
    scale = -scale;
  eta *= scale;

  // TODO: RIGIDBODY
  real sc = scale;
  #pragma acc parallel loop independent deviceptr(vx,vy,vz)
  for (int i = 0; i < n; ++i) {
    vx[i] *= sc;
    vy[i] *= sc;
    vz[i] *= sc;
  }
}
}
TINKER_NAMESPACE_END
