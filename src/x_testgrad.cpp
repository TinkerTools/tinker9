#include "energy.h"
#include "io_text.h"
#include "tinker_rt.h"

TINKER_NAMESPACE_BEGIN
void x_testgrad(int argc, char** argv)
{
   TINKER_RT(initial)();
   TINKER_RT(getxyz)();
   TINKER_RT(mechanic)();

   int flags = calc::xyz;
   flags += (calc::energy + calc::grad);

   rc_flag = flags;
   initialize();
   energy_potential(rc_flag);

   std::vector<real> gdx(n), gdy(n), gdz(n);
   device_array::copyout(n, gdx.data(), gx);
   device_array::copyout(n, gdy.data(), gy);
   device_array::copyout(n, gdz.data(), gz);

   const char* fmt_eng = "\n Total Potential Energy :{:24.4f} Kcal/mole\n\n";
   print(stdout, fmt_eng, esum);

   const char* fmt = " Anlyt{:>10d}       {:12.4f}{:12.4f}{:12.4f}{:14.4f}\n";
   auto do_print = [](int i, int n, int top_m) {
      if (n <= 2 * top_m)
         return true;
      else if (i < top_m)
         return true;
      else if (i >= n - top_m)
         return true;
      else
         return false;
   };
   int print_top_n = 15;
   for (int i = 0; i < n; ++i) {
      if (!do_print(i, n, print_top_n))
         continue;

      real x1 = gdx[i];
      real y1 = gdy[i];
      real z1 = gdz[i];
      real norm = x1 * x1 + y1 * y1 + z1 * z1;
      norm = REAL_SQRT(norm);
      print(stdout, fmt, i + 1, x1, y1, z1, norm);
   }

   finish();
   TINKER_RT(final)();
}
TINKER_NAMESPACE_END
