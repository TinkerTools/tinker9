#include "ff/energy.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/routines.h>

#include "tinker9.h"

#define TINKER_TESTGRAD_VIRIAL 0

namespace tinker {
void xTestgrad(int, char**)
{
   initial();
   tinker_f_getxyz();
   tinker_f_mechanic();
   mechanic2();

   auto out = stdout;
   int digits = inform::digits;

   int flags = (calc::xyz + calc::mass);
   flags += (calc::energy + calc::grad);
#if TINKER_TESTGRAD_VIRIAL
   flags += calc::virial;
#endif

   rc_flag = flags;
   initialize();
   energy(rc_flag);

   std::vector<double> gdx(n), gdy(n), gdz(n);
   copyGradient(calc::grad, gdx.data(), gdy.data(), gdz.data());

   const int len_e = 20 + digits;
   const char* fmt_e = "\n Total Potential Energy :%1$*2$.*3$f Kcal/mole\n\n";
   print(out, fmt_e, esum, len_e, digits);

#if TINKER_TESTGRAD_VIRIAL
   const char* fmt_v = " %-36s%12.3f %12.3f %12.3f\n";
   print(out, fmt_v, "Internal Virial Tensor :", vir[0], vir[1], vir[2]);
   print(out, fmt_v, "", vir[3], vir[4], vir[5]);
   print(out, fmt_v, "", vir[6], vir[7], vir[8]);
#endif

   std::string fmt_t;
   std::string fmt;
   if (digits == 8) {
      fmt_t = "\n  Type    Atom %1$8s "
              "dE/dX %1$9s dE/dY %1$9s dE/dZ %1$9s Norm\n";
      fmt = "\n Anlyt%8d %16.8f%16.8f%16.8f%16.8f";
   } else if (digits == 6) {
      fmt_t = "\n  Type      Atom %1$9s "
              "dE/dX %1$7s dE/dY %1$7s dE/dZ %1$9s Norm\n";
      fmt = "\n Anlyt%10d   %14.6f%14.6f%14.6f  %14.6f";
   } else {
      fmt_t = "\n  Type      Atom %1$12s "
              "dE/dX %1$5s dE/dY %1$5s dE/dZ %1$8s Norm\n";
      fmt = "\n Anlyt%10d       %12.4f%12.4f%12.4f  %12.4f";
   }
   print(out, fmt_t, "");
   // auto do_print = [](int i, int n, int top_m) {
   //    if (n <= 2 * top_m)
   //       return true;
   //    else if (i < top_m)
   //       return true;
   //    else if (i >= n - top_m)
   //       return true;
   //    else
   //       return false;
   // };
   // int print_top_n = 15;
   for (int i = 0; i < n; ++i) {
      // if (not do_print(i, n, print_top_n))
      //    continue;

      real x1 = gdx[i];
      real y1 = gdy[i];
      real z1 = gdz[i];
      real norm = x1 * x1 + y1 * y1 + z1 * z1;
      norm = std::sqrt(norm);
      print(out, fmt, i + 1, x1, y1, z1, norm);
   }

   double totnorm = 0;
   for (int i = 0; i < n; ++i) {
      totnorm += gdx[i] * gdx[i] + gdy[i] * gdy[i] + gdz[i] * gdz[i];
   }
   totnorm = std::sqrt(totnorm);
   double rms = totnorm / std::sqrt(n);
   print(out, "\n\n Total Gradient Norm and RMS Gradient per Atom :\n");
   const char* fmt3 = "\n Anlyt      %1$-30s%2$*3$.*4$f\n";
   const int len3 = 13 + digits;
   print(out, fmt3, "Total Gradient Norm Value", totnorm, len3, digits);
   print(out, fmt3, "RMS Gradient over All Atoms", rms, len3, digits);

   finish();
   tinker_f_final();
}
}
