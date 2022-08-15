#include "ff/amoeba/induce.h"
#include "ff/aplus/induce.h"
#include "ff/atom.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "tool/externfunc.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/polpot.hh>
#include <tinker/detail/units.hh>

namespace tinker {
TINKER_FVOID2(acc1, cu1, diagPrecond2, const real (*)[3], real (*)[3]);
void diagPrecond2(const real (*rsd)[3], real (*zrsd)[3])
{
   TINKER_FCALL2(acc1, cu1, diagPrecond2, rsd, zrsd);
}

void sparsePrecondBuild2() {}

TINKER_FVOID2(acc1, cu1, sparsePrecondApply2, const real (*)[3], real (*)[3]);
void sparsePrecondApply2(const real (*rsd)[3], real (*zrsd)[3])
{
   TINKER_FCALL2(acc1, cu1, sparsePrecondApply2, rsd, zrsd);
}
}

namespace tinker {
TINKER_FVOID2(acc1, cu1, induceMutualPcg2, real (*)[3]);
static void induceMutualPcg2(real (*uind)[3])
{
   TINKER_FCALL2(acc1, cu1, induceMutualPcg2, uind);
}

TINKER_FVOID2(acc1, cu1, induceMutualPcg3, real (*)[3]);
static void induceMutualPcg3(real (*uind)[3])
{
   TINKER_FCALL2(acc1, cu1, induceMutualPcg3, uind);
}

TINKER_FVOID2(acc1, cu1, induceMutualPcg4, real (*)[3]);
static void induceMutualPcg4(real (*uind)[3])
{
   TINKER_FCALL2(acc1, cu1, induceMutualPcg4, uind);
}

void induce2(real (*ud)[3])
{
   if (polpot::use_tholed) {
      induceMutualPcg3(ud);
      ulspredSave(ud, nullptr);
   } else {
      if (polpot::use_expol)
         induceMutualPcg4(ud);
      else
         induceMutualPcg2(ud);
      ulspredSave(ud, nullptr);
   }
   inducePrint(ud);
}
}
