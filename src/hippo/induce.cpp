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
TINKER_F2VOID(cu, 0, acc, 1, diagPrecond2, const real (*)[3], real (*)[3]);
void diagPrecond2(const real (*rsd)[3], real (*zrsd)[3])
{
   TINKER_F2CALL(cu, 0, acc, 1, diagPrecond2, rsd, zrsd);
}

void sparsePrecondBuild2() {}

TINKER_F2VOID(cu, 1, acc, 1, sparsePrecondApply2, const real (*)[3], real (*)[3]);
void sparsePrecondApply2(const real (*rsd)[3], real (*zrsd)[3])
{
   TINKER_F2CALL(cu, 1, acc, 1, sparsePrecondApply2, rsd, zrsd);
}

TINKER_F2VOID(cu, 0, acc, 1, ulspredSave2, const real (*)[3]);
void ulspredSave2(const real (*uind)[3])
{
   TINKER_F2CALL(cu, 0, acc, 1, ulspredSave2, uind);
}

TINKER_F2VOID(cu, 0, acc, 1, ulspredSum2, real (*)[3]);
void ulspredSum2(real (*uind)[3])
{
   TINKER_F2CALL(cu, 0, acc, 1, ulspredSum2, uind);
}
}

namespace tinker {
TINKER_F2VOID(cu, 1, acc, 1, induceMutualPcg2, real (*)[3]);
static void induceMutualPcg2(real (*uind)[3])
{
   TINKER_F2CALL(cu, 1, acc, 1, induceMutualPcg2, uind);
}

TINKER_F2VOID(cu, 1, acc, 1, induceMutualPcg3, real (*)[3]);
static void induceMutualPcg3(real (*uind)[3])
{
   TINKER_F2CALL(cu, 1, acc, 1, induceMutualPcg3, uind);
}

void induce2(real (*ud)[3])
{
   if (polpot::use_dirdamp) {
      induceMutualPcg3(ud);
      ulspredSave3(ud);
   } else {
      induceMutualPcg2(ud);
      ulspredSave2(ud);
   }
   inducePrint(ud);
}
}
