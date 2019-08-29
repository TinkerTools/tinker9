#include "array.h"
#include "energy.h"
#include "io_text.h"
#include "timer.h"
#include "tinker_rt.h"

TINKER_NAMESPACE_BEGIN
static bool timing = true;
// static bool timing = false;

void x_testgrad(int argc, char** argv) {
  if (timing)
    stopwatch_start();

  TINKER_RT(initial)();
  TINKER_RT(getxyz)();
  TINKER_RT(mechanic)();
  if (timing)
    stopwatch_lap("initialized libtinker");

  int flags = calc::xyz;
  flags += (calc::energy + calc::grad);

  rc_flag = flags;
  initialize();
  if (timing)
    stopwatch_lap("initialized libtinkergpu");

  energy_potential(rc_flag & calc::grad);
  if (timing)
    stopwatch_lap("gradient evaluation");

  std::vector<real> gdx(n), gdy(n), gdz(n);
  copyout_array(gdx.data(), gx, n);
  copyout_array(gdy.data(), gy, n);
  copyout_array(gdz.data(), gz, n);
  if (timing)
    stopwatch_lap("gradient copied out");

  const char* fmt = " Anlyt{:>10d}       {:12.4f}{:12.4f}{:12.4f}{:14.4f}\n";
  auto do_print = [](int i, int n) {
    if (n <= 10)
      return true;
    else if (i < 5)
      return true;
    else if (i >= n - 5)
      return true;
    else
      return false;
  };
  for (int i = 0; i < n; ++i) {
    if (!do_print(i, n))
      continue;

    real x1 = gdx[i];
    real y1 = gdy[i];
    real z1 = gdz[i];
    real norm = x1 * x1 + y1 * y1 + z1 * z1;
    norm = REAL_SQRT(norm);
    print(stdout, fmt, i + 1, x1, y1, z1, norm);
  }

  if (timing) {
    stopwatch_stop();
    stopwatch_reset();
  }
}
TINKER_NAMESPACE_END
