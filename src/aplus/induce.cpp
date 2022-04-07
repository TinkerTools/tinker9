#include "ff/amoeba/induce.h"
#include "ff/hippo/induce.h"
#include "ff/nblist.h"

namespace tinker {
void diagPrecond3(const real (*rsd)[3], real (*zrsd)[3])
{
   diagPrecond2(rsd, zrsd);
}

void sparsePrecondBuild3() {}

extern void sparsePrecondApply3_acc(const real (*)[3], real (*)[3]);
extern void sparsePrecondApply3_cu(const real (*)[3], real (*)[3]);
void sparsePrecondApply3(const real (*rsd)[3], real (*zrsd)[3])
{
#if TINKER_CUDART
   if (ulistVersion() & Nbl::SPATIAL)
      sparsePrecondApply3_cu(rsd, zrsd);
   else
#endif
      sparsePrecondApply3_acc(rsd, zrsd);
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

namespace tinker {
extern void induceMutualPcg3_acc(real (*uind)[3]);
extern void induceMutualPcg3_cu(real (*uind)[3]);
static void induceMutualPcg3(real (*uind)[3])
{
#if TINKER_CUDART
   if (pltfm_config & Platform::CUDA)
      induceMutualPcg3_cu(uind);
   else
#endif
      induceMutualPcg3_acc(uind);
}

void induce3(real (*uind)[3])
{
   induceMutualPcg3(uind);
   ulspredSave3(uind);
   inducePrint(uind);
}
}
