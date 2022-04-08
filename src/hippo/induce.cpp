#include "ff/amoeba/induce.h"
#include "ff/aplus/induce.h"
#include "ff/atom.h"
#include "ff/nblist.h"
#include "ff/potent.h"
#include "tool/ioprint.h"
#include <tinker/detail/inform.hh>
#include <tinker/detail/polar.hh>
#include <tinker/detail/polpot.hh>
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
