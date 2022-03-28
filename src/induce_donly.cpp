#include "ff/hippo/inducechgpen.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "md/inc.h"
#include "tool/io.h"
#include <cmath>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/potent.hh>
#include <tinker/detail/units.hh>
#include <tinker/detail/uprior.hh>

namespace tinker {
void sparse_precond_apply2_acc(const real (*)[3], real (*)[3]);
void sparse_precond_apply2_cu(const real (*)[3], real (*)[3]);
void sparsePrecondBuild2() {}

void sparsePrecondApply2(const real (*rsd)[3], real (*zrsd)[3])
{
#if TINKER_CUDART
   if (ulistVersion() & Nbl::SPATIAL)
      sparse_precond_apply2_cu(rsd, zrsd);
   else
#endif
      sparse_precond_apply2_acc(rsd, zrsd);
}

void ulspred_save2_acc(const real (*)[3]);
void ulspred_sum2_acc(real (*)[3]);
void ulspredSave2(const real (*uind)[3])
{
   ulspred_save2_acc(uind);
}

void ulspredSum2(real (*uind)[3])
{
   ulspred_sum2_acc(uind);
}

void induce_mutual_pcg2(real (*uind)[3]);
void induce_mutual_pcg2_acc(real (*uind)[3]);
void induce_mutual_pcg2_cu(real (*uind)[3]);

void induce2(real (*ud)[3])
{
   induce_mutual_pcg2(ud);
   ulspredSave2(ud);

   if (inform::debug && usePotent(Potent::POLAR)) {
      std::vector<double> uindbuf;
      uindbuf.resize(3 * n);
      darray::copyout(g::q0, n, uindbuf.data(), ud);
      wait_for(g::q0);
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
