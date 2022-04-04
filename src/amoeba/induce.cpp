#include "ff/atom.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/units.hh>

namespace tinker {
extern void diagPrecond_acc(const real (*rsd)[3], const real (*rsdp)[3], //
   real (*zrsd)[3], real (*zrsdp)[3]);
void diagPrecond(const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3], real (*zrsdp)[3])
{
   diagPrecond_acc(rsd, rsdp, zrsd, zrsdp);
}

void sparsePrecondBuild() {}

extern void sparsePrecondApply_acc(const real (*)[3], const real (*)[3], real (*)[3], real (*)[3]);
extern void sparsePrecondApply_cu(const real (*)[3], const real (*)[3], real (*)[3], real (*)[3]);
void sparsePrecondApply(
   const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3], real (*zrsdp)[3])
{
#if TINKER_CUDART
   if (ulistVersion() & Nbl::SPATIAL)
      sparsePrecondApply_cu(rsd, rsdp, zrsd, zrsdp);
   else
#endif
      sparsePrecondApply_acc(rsd, rsdp, zrsd, zrsdp);
}

extern void ulspredSave_acc(const real (*)[3], const real (*)[3]);
void ulspredSave(const real (*uind)[3], const real (*uinp)[3])
{
   ulspredSave_acc(uind, uinp);
}

extern void ulspredSum_acc(real (*)[3], real (*)[3]);
void ulspredSum(real (*uind)[3], real (*uinp)[3])
{
   ulspredSum_acc(uind, uinp);
}
}

namespace tinker {
extern void induceMutualPcg1_acc(real (*uind)[3], real (*uinp)[3]);
extern void induceMutualPcg1_cu(real (*uind)[3], real (*uinp)[3]);
static void induceMutualPcg1(real (*uind)[3], real (*uinp)[3])
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      induceMutualPcg1_cu(uind, uinp);
   else
#endif
      induceMutualPcg1_acc(uind, uinp);
}

void induce(real (*ud)[3], real (*up)[3])
{
   induceMutualPcg1(ud, up);
   ulspredSave(ud, up);

   if (inform::debug and usePotent(Potent::POLAR)) {
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
}
