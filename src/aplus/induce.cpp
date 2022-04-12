#include "ff/amoeba/induce.h"
#include "ff/hippo/induce.h"
#include "ff/nblist.h"
#include "tool/externfunc.h"

namespace tinker {
void diagPrecond3(const real (*rsd)[3], real (*zrsd)[3])
{
   diagPrecond2(rsd, zrsd);
}

void sparsePrecondBuild3() {}

TINKER_FVOID2(cu, 1, acc, 1, sparsePrecondApply3, const real (*)[3], real (*)[3]);
void sparsePrecondApply3(const real (*rsd)[3], real (*zrsd)[3])
{
   TINKER_FCALL2(cu, 1, acc, 1, sparsePrecondApply3, rsd, zrsd);
}

void ulspredSave3(const real (*uind)[3])
{
   ulspredSave2(uind);
}

void ulspredSum3(real (*uind)[3])
{
   ulspredSum2(uind);
}
}
