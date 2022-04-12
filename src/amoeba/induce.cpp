#include "ff/atom.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "tool/externfunc.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/units.hh>

namespace tinker {
TINKER_FVOID2(
   cu, 0, acc, 1, diagPrecond, const real (*)[3], const real (*)[3], real (*)[3], real (*)[3]);
void diagPrecond(const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3], real (*zrsdp)[3])
{
   TINKER_FCALL2(cu, 0, acc, 1, diagPrecond, rsd, rsdp, zrsd, zrsdp);
}

void sparsePrecondBuild() {}

TINKER_FVOID2(cu, 1, acc, 1, sparsePrecondApply, const real (*)[3], const real (*)[3], real (*)[3],
   real (*)[3]);
void sparsePrecondApply(
   const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3], real (*zrsdp)[3])
{
   TINKER_FCALL2(cu, 1, acc, 1, sparsePrecondApply, rsd, rsdp, zrsd, zrsdp);
}

TINKER_FVOID2(cu, 0, acc, 1, ulspredSave, const real (*)[3], const real (*)[3]);
void ulspredSave(const real (*uind)[3], const real (*uinp)[3])
{
   TINKER_FCALL2(cu, 0, acc, 1, ulspredSave, uind, uinp);
}

TINKER_FVOID2(cu, 0, acc, 1, ulspredSum, real (*)[3], real (*)[3]);
void ulspredSum(real (*uind)[3], real (*uinp)[3])
{
   TINKER_FCALL2(cu, 0, acc, 1, ulspredSum, uind, uinp);
}
}

namespace tinker {
TINKER_FVOID2(cu, 1, acc, 1, induceMutualPcg1, real (*)[3], real (*)[3]);
static void induceMutualPcg1(real (*uind)[3], real (*uinp)[3])
{
   TINKER_FCALL2(cu, 1, acc, 1, induceMutualPcg1, uind, uinp);
}

void inducePrint(const real (*ud)[3])
{
   if (inform::debug and use(Potent::POLAR)) {
      std::vector<double> uindbuf;
      uindbuf.resize(3 * n);
      darray::copyout(g::q0, n, uindbuf.data(), ud);
      waitFor(g::q0);
      bool header = true;
      for (int i = 0; i < n; ++i) {
         if (polar::polarity[i] != 0) {
            if (header) {
               header = false;
               print(stdout, "\n Induced Dipole Moments (Debye) :\n");
               print(stdout, "\n    Atom %1$13s X %1$10s Y %1$10s Z %1$9s Total\n\n", "");
            }
            double u1 = uindbuf[3 * i];
            double u2 = uindbuf[3 * i + 1];
            double u3 = uindbuf[3 * i + 2];
            double unorm = std::sqrt(u1 * u1 + u2 * u2 + u3 * u3);
            u1 *= units::debye;
            u2 *= units::debye;
            u3 *= units::debye;
            unorm *= units::debye;
            print(stdout, "%8d     %13.4f%13.4f%13.4f %13.4f\n", i + 1, u1, u2, u3, unorm);
         }
      }
   }
}

void induce(real (*ud)[3], real (*up)[3])
{
   induceMutualPcg1(ud, up);
   ulspredSave(ud, up);
   inducePrint(ud);
}
}
