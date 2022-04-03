#include "ff/amoeba/induce.h"
#include "ff/nblist.h"

namespace tinker {
void induce_mutual_pcg1_acc(real (*uind)[3], real (*uinp)[3]);
void induce_mutual_pcg1_cu(real (*uind)[3], real (*uinp)[3]);

void sparse_precond_apply_acc(const real (*)[3], const real (*)[3], real (*)[3], real (*)[3]);
void sparse_precond_apply_cu(const real (*)[3], const real (*)[3], real (*)[3], real (*)[3]);

void ulspred_save_acc(const real (*)[3], const real (*)[3]);
void ulspred_sum_acc(real (*)[3], real (*)[3]);
}

namespace tinker {
void sparsePrecondBuild() {}

void sparsePrecondApply(
   const real (*rsd)[3], const real (*rsdp)[3], real (*zrsd)[3], real (*zrsdp)[3])
{
#if TINKER_CUDART
   if (ulistVersion() & Nbl::SPATIAL)
      sparse_precond_apply_cu(rsd, rsdp, zrsd, zrsdp);
   else
#endif
      sparse_precond_apply_acc(rsd, rsdp, zrsd, zrsdp);
}

void ulspredSave(const real (*uind)[3], const real (*uinp)[3])
{
   ulspred_save_acc(uind, uinp);
}

void ulspredSum(real (*uind)[3], real (*uinp)[3])
{
   ulspred_sum_acc(uind, uinp);
}
}

#include "ff/amoeba/elecamoeba.h"
#include "ff/amoeba/empole.h"
#include "ff/amoeba/epolar.h"
#include "ff/amoeba/induce.h"
#include "ff/elec.h"
#include "ff/energy.h"
#include "ff/nblist.h"
#include "ff/pme.h"
#include "ff/potent.h"
#include "math/zero.h"
#include "tool/io.h"
#include <cmath>
#include <map>
#include <tinker/detail/couple.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mplpot.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/polgrp.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/sizes.hh>
#include <tinker/detail/units.hh>
#include <tinker/detail/uprior.hh>
namespace tinker {
extern void induce_mutual_pcg1(real (*uind)[3], real (*uinp)[3]);
void induce(real (*ud)[3], real (*up)[3])
{
   induce_mutual_pcg1(ud, up);
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
