#include "ff/atom.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/units.hh>

namespace tinker {
extern void diagPrecond2_acc(const real (*rsd)[3], real (*zrsd)[3]);
void diagPrecond2(const real (*rsd)[3], real (*zrsd)[3])
{
   diagPrecond2_acc(rsd, zrsd);
}

void sparsePrecondBuild2() {}

extern void sparsePrecondApply2_acc(const real (*)[3], real (*)[3]);
extern void sparsePrecondApply2_cu(const real (*)[3], real (*)[3]);
void sparsePrecondApply2(const real (*rsd)[3], real (*zrsd)[3])
{
#if TINKER_CUDART
   if (ulistVersion() & Nbl::SPATIAL)
      sparsePrecondApply2_cu(rsd, zrsd);
   else
#endif
      sparsePrecondApply2_acc(rsd, zrsd);
}

extern void ulspredSave2_acc(const real (*)[3]);
void ulspredSave2(const real (*uind)[3])
{
   ulspredSave2_acc(uind);
}

extern void ulspredSum2_acc(real (*)[3]);
void ulspredSum2(real (*uind)[3])
{
   ulspredSum2_acc(uind);
}
}

namespace tinker {
extern void induceMutualPcg2_acc(real (*uind)[3]);
extern void induceMutualPcg2_cu(real (*uind)[3]);
static void induceMutualPcg2(real (*uind)[3])
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      induceMutualPcg2_cu(uind);
   else
#endif
      induceMutualPcg2_acc(uind);
}

void induce2(real (*ud)[3])
{
   induceMutualPcg2(ud);
   ulspredSave2(ud);

   if (inform::debug && use(Potent::POLAR)) {
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
